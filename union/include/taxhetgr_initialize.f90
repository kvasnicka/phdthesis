!Set the number of threads to be used by MKL (when we use all of physical cores to distribute workload
!evenly, it doesn't make sense to try to set this to a higher value than 1). This setting does not
!matter if the compiler option is -mkl=sequential.
call mkl_set_num_threads(use_threads_mkl)

!CAF variables
i_img = this_image() !index of current image
n_img = num_images() !number of images


!Tracking runtime
initial_time = time_real_acc()

if(I_par/=2) then
    write(*,*) 'Error: At this stage many of the parts of programme assume I=2.'
    error stop
end if

!I/O handled by one image only.
if(i_img == spec_img) then
    !Initialize output and read parameters from file
    parameter_file = get_parfilename() !Get the name of the parameter file (from command line argument)

    call initialize_output(folder_name,parameter_file) !create folder for saving results, etc.

    call read_parameters(parameter_file,pars) !read parameters

    call check_parameters(pars,folder_name)

    call double_output(folder_name,'****************************************************')
    call double_output(folder_name,'Optimal policy in a union')
    call double_output(folder_name,'****************************************************')

!The following is useful only on HPC in slurm output (hence not saved in logfile)
write(*,*) 'Parameter file: ',parameter_file
end if

!If IID_shocks, then the number of grid points for shock realization will be equal to one.
if(pars%IID_shocks) then
    pars%M_grid = 1
else
    pars%M_grid = M_par
end if

!Compute bounds for maximization over choice vector x (as defined in x_guess).
!(This only uses resource constraint and tie endowment - of course extreme
!choices in this regard will never be optimal, but it's better not to rule them out
!at this stage)
!First get consumption matrix (max feasible choice),

call max_feas_cons(pars%l_max,pars%mass,pars%theta_0,pars%G,pars%max_cons_mat)
!First M elements of x_bu are the labour endowment
pars%x_bu(1:M_par) = reshape(pars%l_max,[M_par])
!The remainder - get it using subroutine  c_inc_mtv(c_inc_mat,c_inc_vec)
call c_inc_mtv(pars%max_cons_mat,pars%x_bu(1+M_par:3*M_par -1))


!retrieve the parameters at all other images
sync all
if(i_img /= spec_img) pars = pars[spec_img]

!Get the Unix time at which VFI algorithm should stop (strictly speaking need it at one
!image only but it doesn't matter and it could lead to bugs later if we tried to use it at some
!other image).
runtime_limit_unix = initial_time + pars%runtime_limit * 3600.0_dp
!spend at most 10% of runtime in CTS algorithm
runtime_limit_unix_CTS = initial_time + (pars%runtime_limit * 3600.0_dp)/10.0_dp
!and at most 80% of the time in Howard (we want to converge smoothly and avoid overshooting)
runtime_limit_unix_How = initial_time + (pars%runtime_limit * 3600.0_dp)*0.8_dp

!Local copy of the parameters (this is a tiny object, no big deal), we create it
!only so we can associate a pointer with its elements, and then access it withing various
!subroutines
pars_img = pars
pars_pntr => pars_img

!If pars%load_initial_guess == .true., then we don't want to use the CTS algorithm
!to obtain an initial guess from a coarser grid problem
if(pars%load_initial_guess == .true. .and. pars%CTS_split_times > 0) then
    pars%CTS_split_times = 0
    if(i_img==spec_img) then
        call double_output(folder_name,'Loading initial guess from file. CTS algorithm not used.')
    end if
end if

!If sim_only == .true., then we want to load Value function from a file
!and proceed with solving the first-period problem and simulation.
!in this case, set CTS_split_times = 0 and set load_initial_guess = .true..
!doing it this way saves us time because we don't have to rewrite code
!for generating grids, etc.
if(pars%sim_only == .true.) then
    pars%CTS_split_times = 0
    pars%load_initial_guess = .true.
    if(i_img==spec_img) then
        call double_output(folder_name,'Simulating only, loading value function from folder input_folder.')
    end if
end if

if(pars%sim_interp == .true. .and. pars%sim_only == .false.) then
    if(this_image() == 1) then
    call double_output(folder_name,'The interpolation of policy function is supported only when sim_only == .true.')
    call double_output(folder_name,'Setting sim_interp = .false.')
    end if
    pars%sim_interp = .false.
    !If we want to add this support then we need to do quite a bit of work in terms of checking what is allocated and
    !what is not, and pass the policy function between images. But it is pretty much useless because we want to
    !run large scale simulations separately from the value function iteration - therefore for now, I implemented
    !this option only for loading the value function from file.
end if

!If we are using Howard's acceleration algorithm, set VFI_max_iter to a really large
!number so it is never reached, and we thus use only the time limit!
if(pars%crit_How > 0) then
    if(this_image()==1) then
    call double_output(folder_name,'Howard algorithm used. VFI_max_iter will not be used.')
    end if
    pars%VFI_max_iter = 1000000
end if

allocate(shep_inc(pars%shep_n**4,4)) !save 4 indices (not M)
!Fill shep_inc with all the possible values of increments
do shep_ind = 1,pars%shep_n**4
    call ind_FtU(shep_ind,[pars%shep_n,pars%shep_n,pars%shep_n,pars%shep_n],shep_inc(shep_ind,:))
end do
!subtract 1 so the increments range from 0 and we don't have to subtract 1 every time we use them
!deep in the program.
shep_inc = shep_inc - 1
shep_inc_pntr => shep_inc

!Also save the parameters
if(this_image() == 1) then
    call save_par(pars,folder_name)
end if

N_sim_done = 0

stop_now = .false.

sync all !sync before starting solutions so that all images have parameters.
