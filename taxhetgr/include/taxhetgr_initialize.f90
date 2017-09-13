!Set the number of threads to be used by MKL (when we use all of physical cores to distribute workload
!evenly, it doesn't make sense to try to set this to a higher value than 1). This setting does not
!matter if the compiler option is -mkl=sequential.
call mkl_set_num_threads(use_threads_mkl)

!CAF variables
i_img = this_image() !index of current image
n_img = num_images() !number of images

!(EFE) move definition of spec_img from parameters_hard here and pass it
!to subroutines where it is needed. spec_img should depend on the number
!of images (at least potentially). A good choice might be
!spec_img = num_images.

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
    call double_output(folder_name,'Taxation with Heterogeneous Growth (TaxHetGr)')
    call double_output(folder_name,'****************************************************')

!The following is useful only on HPC in slurm output (hence not saved in logfile)
write(*,*) 'Parameter file: ',parameter_file
end if

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
!only so we can associate a pointer with its elements.
pars_img = pars
!associate pointers (we don't want to point to the whole object pars because in the future
!it may contain allocatable arrays)
S_pntr => pars_img%S
P_pntr => pars_img%P
A_par_pntr => pars_img%A_par
B_par_pntr => pars_img%B_par
Gamma_par_pntr => pars_img%Gamma_par
mass_pntr => pars_img%mass

!This defeats the purpose of the things above and is temporary. Of course it doesn't
!make sense to point only to certain elements of pars because in the future
!it might contain allocatable arrays, and then point to the whole thing anyway.
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

!If we are using Howard's acceleration algorithm, set VFI_max_iter to a really large
!number so it is never reached, and we thus use only the time limit!
if(pars%crit_How > 0) then
    if(this_image()==1) then
    call double_output(folder_name,'Howard algorithm used. VFI_max_iter will not be used.')
    end if
    pars%VFI_max_iter = 1000000
end if

sync all !sync before starting solutions so that all images have parameters.
