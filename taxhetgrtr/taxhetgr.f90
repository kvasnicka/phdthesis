program taxhetgr
!This is the main file of the program taxhetgrtr (taxation under heterogeneous growth with transfers)



use mod_types !types for integer/real variables
use mod_IO !subroutines used for input/output
use mod_par !definition of type par and subroutine for loading this from parameter file
use mod_tools !various tools (tracking runtime, etc.)
use mod_parameters_hard !parameters which are hard-coded and not loaded from parameter file
use mod_work_dist !tools for distributing workload between images (processes)
use mod_taxhetgr !types, subroutines, functions useful in this project only
use mod_utility

use nag_library

use ifport !so we can use system()

implicit none

!___________________________________________________________________
!Variable declaration: later on either put this into an include statement or a module
character(256) :: parameter_file = ''

character(256) :: folder_name = ''
character(256) :: elapsed_time_string
type(par), codimension[*] :: pars !Type which contains all parameters, defined in module mod_par
!this can either be a coarray, or reading of data from parameter file could be done
type(par), target :: pars_img !Local copy of parameters, used because we need to associated its elements
!to a pointer and we can't give the TARGET attribute to a coarray.

real(dp) :: initial_time !Initial unix time including miliseconds.

!CAF control variables
integer :: i_img = 0
integer :: n_img = 0
integer :: V_num_el = 0 !number of elements of the value function.
integer :: min_ind_img = 0
integer :: max_ind_img = 0

integer :: min_ind_sim_img = 0
integer :: max_ind_sim_img = 0
integer :: N_sim_loc !local version of pars%N_sim (will be equal to 1 on all but one image - used to save memory
!when allocating - just a technicality due to use of coarrays)

!
!Value function
real(dp), allocatable, dimension(:,:,:,:), target :: V_old !Allocated all the time at all images, 4 states if I=2.
real(dp), allocatable, dimension(:) :: V_stub !Piece of value function unfolded into row vector, allocated on all images

type(shared_data_V_type), codimension[*] :: shared_data_V
!!allocated on special image (index spec_img) only, used for passing results
!!(more precisely, this is allocated on all images, but its components are allocated to size 1
!!on all images but one). Also note that reallocating components of this doesn't lead to implicit
!!synchronization as reallocating whole coarrays does -> need to sync all manually!!!
!
!Some VFI algorithm iteration-invariant data saved on every image:
integer, allocatable, dimension(:,:) :: ind_unfold_all !contains all of the 'unfolded' indices for every image,
!where unfolded index refers to the value function in the multi-dimensional as opposed the vector form.
!This is just to save time when accessing data.


logical, codimension[*] :: VFI_stop_now !This value will be determined on image spec_img only
!all images will stop if it is .true.
!
!CTS algorithm variables
integer :: CTS_ind = 0
!
!Initial guess and its sharing (later on could be moved to a module)
type(c_guess_all), allocatable, dimension(:) :: c_guess_img !initial consumption guess for all grid points at which
!the image works (division b/w images the same as division of value function)

!Index of optimal value of what_binds at all gridpoints on an images.
integer, dimension(:), allocatable :: optim_wb_img

!The policy function at a grid point. C_pol_stub(:,:,j) contains
!consumption matrix which is optimal at grid point with index j.
!Also save the old function so we can check convergence in terms of policy function...
real(dp), allocatable, dimension(:,:,:) :: C_pol_stub,C_pol_stub_old
real(dp), dimension(M_par,I_par) :: C_pol_gp !optimal choice at a grid point

real(dp), codimension[*] :: share_feas_img !weighted share of grid points at which feasible solution was found.
real(dp), codimension[*] :: avg_loss_img !average loss at image (loss from violating constraints)
integer :: img_ind !used when we need to cycle over images
integer :: sim_ind !used when we need to cycle over samples in the simulation
!
!!variables for generating output (useful when we need to save formatted output).
character(256) :: output_string, tmp_string1,tmp_string2,tmp_string3
!
!Numbers of grid points (these are not the same as those in the parameter file because they can be
!changed by CTS algorithm)
integer :: N_a = 0 !number of grid points for asset holdings
integer :: N_rho = 0 !number of grid points for ratio of marginal utilities
integer :: N_T = 0 !number of grid points for 'time' governing deterministic trends
!(number of grid points per shock process is given in parameters_hard: M_par)


!!Temporary variables: used for loading the value function from file and in between
!!iterations of CTS algorithm.
integer :: N_a_tmp,N_rho_tmp,N_t_tmp
real(dp), dimension(:), allocatable :: V_vect_tmp !temporary value function in vector form (used when loading
!value function from file and in CTS algorithm).
!
!real(dp) :: share_feas,avg_loss
!
integer :: V_ind !VFI algorithm index

type(grids_type), target :: grids, grids_tmp !grids and temporary grids (the latter used in CTS algorithm and
!when loading a value function from file.

!Variables used in search for feasible solution.
integer :: num_fc
real(dp), dimension(:,:,:), allocatable :: LFFC_choices, LFFC_choices_all
!using 2 arrays for essentially the same data because we want to save some memory
!by copying the result of (ex ante) unknown size into a smaller array

!!index for cycling over grid points
integer :: gp_ind

!!index for cycling over images
!integer :: img_ind
!
logical :: skip_max = .false. !If this is true, the maximization step will be skipped
!in the current iteration of VFi algorithm (Howard)
real(dp) :: V_new_gp !new value of V at a grid point
logical :: glob_opt

integer :: rng_seed = 0 !some value for initializing RNG

!integer(4) :: sys_call_result
!
!runtime_limit_unix: unix time at which VFI should stop
real(dp) :: runtime_limit_unix
real(dp) :: runtime_limit_unix_CTS !runtime limit for CTS algorithm
!real(dp) :: time_unix

integer :: s_init_ind_a !actual initial shock

!Pointer to value function. This will be used to make the old value function V_old
!accessible to the subroutine used in NAG maximization.
real(dp), dimension(:,:,:,:), pointer :: V_old_pntr !4 states when I=2
!Pointers for passing grids. Unfortunately, the type grids_type contains
!allocatable arrays and pointers to the whole type could not be part of a
!common block. Therefore I define pointers for each component of the grid.

real(dp), dimension(:), pointer :: a1_gr_pntr,rho_gr_pntr,t_gr_pntr
real(dp), dimension(:,:), pointer :: S_pntr,P_pntr !points to shock space
!and transition matrix.
real(dp), pointer :: A_par_pntr,B_par_pntr,Gamma_par_pntr
real(dp), dimension(:,:), pointer :: mass_pntr
type(par), pointer :: pars_pntr !This will be viable only as long as pars contain
!no allocatable arrays. Right now it is easier to do this than reweite get_c_im
!and the function which it uses but when I generalize to I>2 and if I want the code
!to work for arbitrary M w/o recompilation, I will need to change this.
!!___________________________________________________________________________
!
!!Maximum (Average) Relative Difference in policy function between
!!iterations. The average is also weighted by the number of grid points.
real(dp), codimension[*] :: C_pol_MARD_img, C_pol_AARD_img
real(dp), codimension[*] :: C_pol_MARD, C_pol_AARD

        !The following variables are use for randomly drawing initial shock realization
        !if s_init_index = -1 in parameter file. We want to draw the shock with frequency associated
        !with stationary distribution.
        real(dp) :: k,r1,r2,randnum

!!Debug variables (temporary)
!integer :: problem_ind !index of grid point where a problem occurs
!integer, dimension(5) :: problem_ind_unfold !unfolded index of a gp where problem occurs
!!end of tmp debug vars
!real(dp) :: debug_val,debug_val2
!
real(dp) :: loss_gp, loss_all_avg !loss at grid point, average loss at all grid points.
real(dp), codimension[*] :: loss_img_tot !loss average across gp's that an image works on

real(dp) :: V_MARD,V_AARD !maximum and average absolute relative difference.
real(dp), dimension(1) :: V_quant_ARD, quant_nag_input
real(dp), dimension(:), allocatable :: ARD_vector !vector of all absolute relative differences
integer :: ifail !for use in NAg
!_____________________________________________________________________
type(FP_result_type), dimension(M_par) :: FP_result !result of the first-period problem for every shock realization

type(sim_result_type) :: sim_result !simulated series

type(sim_result_type) :: sim_avg !sample average of all simulated series (tbh this could just
!be part of sim_res_all)

integer :: m_ind

!_____________________________________________________________________
!character(MAX_HOSTNAM_LENGTH + 1) :: host_name !need USE IFPORT for this
!!This is used in debugging only to find out which node an image runs on
!
logical, codimension[*] :: use_How !is true if using Howard acceleration algorithm. Will be set to .false.
!once a criterion is less than c_How (see parameter file for definition and discussion)
integer, codimension[*] :: skip_count_How !counter of maximization steps to be skipped
real(dp) :: crit_value_How !value of the criterion
integer :: i_How

!testing:
type(sim_result_all_type) :: sim_res_all !This contains all the simulated series (will be allocated on image 1 only),
!and some other results.


common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
!to some subroutines.
common /pars_pointer/ pars_pntr

common /comm1/ a1_gr_pntr,rho_gr_pntr,t_gr_pntr,A_par_pntr,B_par_pntr,Gamma_par_pntr,S_pntr,P_pntr,mass_pntr
!!___________________________________________________________________________


include 'include/taxhetgr_initialize.f90'

!!______________________________________________________________________________
do CTS_ind = 1,pars%CTS_split_times + 1


!!Howard acceleration algorithm - at every stage of CTS initialize the variable telling the
!!program whether to use acceleration.
if(pars%crit_How /= 0) then
    use_How = .true.
else
    use_How = .false.
    if (i_img == spec_img) call double_output(folder_name,'Using Howard acceleration algorithm')
end if


skip_count_How = 0 !Counter of maximization steps skipped - will be set to some positive value
!whenever a criterion is satisfied and every time a maximization step is skipped this counter
!will decrease by one.
skip_max = .false.


!Compute number of grid points at this iteration of CTS algorithm (at every grid point). If CTS
!algorithm is not used, then the values in parameter file will be returned by the subroutine.
call CTS_num_points(CTS_ind,pars,N_a,N_rho,N_T)
V_num_el = N_a*N_rho*N_T*M_par !number of elements of value function (used often)

!Display information about CTS algorithm and workload distribution.
if(i_img == spec_img .and. .not.(pars%sim_only)) then !Write this even if CTs not used (then CTS_ind = 1 out of 1)
    write(tmp_string1,'(I3)') CTS_ind
    write(tmp_string2,'(I3)') pars%CTS_split_times + 1
    output_string = 'CTS algorithm iteration '//trim(tmp_string1) //' out of '//trim(tmp_string2)
    call double_output(folder_name,trim(output_string))

    write(tmp_string1,'(I6)') N_a
    write(tmp_string2,'(I6)') N_rho
    write(tmp_string3,'(I6)') N_t

    output_string = 'Number of gridpoints: N_a = '//trim(tmp_string1)//', N_rho = '//trim(tmp_string2)// &
    ', N_T = '//trim(tmp_string3)
    call double_output(folder_name,trim(output_string))

    write(tmp_string1,'(I6)') n_img
    write(tmp_string2,'(F10.1)') real(V_num_el)/real(n_img)
    output_string = 'Working on '//trim(tmp_string1)//' images ('//trim(tmp_string2)//' gridpoints per image)'
    call double_output(folder_name,trim(output_string))
    call double_output(folder_name,'____________________________________________________')
end if


!Get the range of indices (in the vector representation of value function)
!on which every image will work.
call get_img_ind(i_img,n_img,V_num_el,min_ind_img,max_ind_img)


!allocate memory to the unfolded indices. (indexed the same way as V_stub, so that
!element with index i contains the unfolded index for i-th term in the vector representation
!of the value function.
allocate(ind_unfold_all(4,min_ind_img:max_ind_img)) !4 state variables -> This already assumes I_par=2



!Get the unfolded indices (indices in the multidimensional value function corresponding to the
!indices in vector form value function)
call get_ind_unfold(min_ind_img,max_ind_img,[N_a,N_rho,N_T,M_par],ind_unfold_all)


!Generate grids:
allocate(grids%A1_gr(N_a),grids%rho_gr(N_rho),grids%T_gr(N_T))
call grids_gen(grids%A1_gr,grids%rho_gr,grids%T_gr,N_a,N_rho,N_T,pars)

!Associate pointers with elements of the grids
a1_gr_pntr => grids%A1_gr
rho_gr_pntr => grids%rho_gr
t_gr_pntr => grids%t_gr



if(pars%sim_only) goto 999 !if we only want to simulate, skip LFFC.


!__________________INITIAL GUESS of consumption___________________________
!Before starting value function iteration, we need to get an initial guess at every grid point.
!To begin with, we look for any feasible choice. In VFI algorithm, the initial guess
!will be updated at every iteration so that we hopefully start close to the optimum and
!save time. In later stages of VFI algorithm, Howard improvement algorithm will be used
!which essentially means that a consumption choice obtained as optimum at an iteration
!of VFI algorithm will be used for several iterations (skipping the maximization step).

!allocate memory!(same indexing convention as in the value function, so indiced
!correspond to elements of the whole value function in vector form, not elements
!of the piece of value function that an image works on.
allocate(c_guess_img(min_ind_img:max_ind_img))

!Also allocate memory to the array which contains indices of optimal what_binds
allocate(optim_wb_img(min_ind_img:max_ind_img))
!initialize to avoid some weird bugs
optim_wb_img = 1


!Get the guesses of consumption at the outset (it is useful to do this this way
!so that it is already passed to subroutine c_guess_img), particularly so that it
!will then be easy to generalize the algorithm for finding feasible grids.
!The guesses here are in form of shares of maximum consumption feasible when we
!only take the resource constraint into account.
allocate(LFFC_choices_all(M_par,I_par,pars%LFFC_max_ind**(I_par*M_par-1)))
call LFFC_prep_choices(pars,LFFC_choices_all,num_fc)

!Now that we know the number of choices found by the algorithm, put them into
!a smaller array and free up some memory
allocate(LFFC_choices(M_par,I_par,num_fc))
LFFC_choices(:,:,1:num_fc) = LFFC_choices_all(:,:,1:num_fc)
deallocate(LFFC_choices_all)


if(i_img == spec_img) then
    call double_output(folder_name,'Looking for feasible choices...')
end if

call get_c_guess_img(c_guess_img,min_ind_img,max_ind_img,ind_unfold_all,pars,grids,&
share_feas_img,avg_loss_img,LFFC_choices,num_fc)

!Sync after every image gets an initial guess of consumption. This is mainly so that we
!can compute the overall feasibility share and display a message about it (every image
!only needs its own copy of c_guess_img, no need to access initial guesses for gridpoints
!which are processed by other images).
sync all

!Display a message about results of LFFC
call LFFC_msg(share_feas_img,avg_loss_img,min_ind_img,max_ind_img,spec_img,V_num_el,folder_name,initial_time)

!After this point we no longer have available the array of all of initial guesses possibilities!
deallocate(LFFC_choices)



!The following section serves to get an initial guess of value function to all images. If CTS
!algorithm is used, and we are at the first iteration, we either load the guess from file
!by one image and send it to all images, or we do the same thing with a constant guess.
!If we are not at the first stage of CTS algorithm, one image computes a value
!function using interpolation from the value function computed on the coarse grid. This
!function is then sent to all images.

!We prepare the data to be sent on one image, then we synchronize and send it to all images.
999 if(i_img==spec_img) then

    !Allocate memory to object where we will save the new value function on one image only
    allocate(shared_data_V%V_tmp(N_a,N_rho,N_t,M_par))


    if(CTS_ind == 1) then
        !Either load an initial guess from input_folder (given in parameter file) or use an initial guess
        !given here (or in an include file if it's something more complicated)
        if(pars%load_initial_guess==.true.) then
            !load the guess from old file on special image

            !Get number of grid points from file
            call load_N(N_a_tmp,N_rho_tmp,N_t_tmp,pars%input_folder,folder_name)
            !allocate memory to the temporary value function and grids and load these from files
            allocate(V_vect_tmp(N_a_tmp*N_rho_tmp*N_t_tmp*M_par))
            allocate(grids_tmp%A1_gr(N_a_tmp),grids_tmp%rho_gr(N_rho_tmp),&
                grids_tmp%t_gr(N_t_tmp))

            call load_V(V_vect_tmp,grids_tmp%a1_gr,grids_tmp%rho_gr, &
                grids_tmp%t_gr,pars%input_folder,folder_name)

            !(EFE) could paralelize the interpolation but there is not much to gain,
            !the code spends very little time here compared to other places.

            !reshaping the vector form value function so that is can be an input for V_from_V
            !This depends on the order of gridpoints in the vector form being correct.

            call V_from_V(reshape(V_vect_tmp,[N_a_tmp,N_rho_tmp,N_t_tmp,M_par]),&
                grids_tmp,grids,shared_data_V%V_tmp,folder_name)

            !Now the interpolated value function is saved in shared_data_V%V_tmp, which can
            !be retrieved by all images.

            !We already have the new value function, do not need to keep the old value
            !function or the old grid
            deallocate(V_vect_tmp)
            deallocate(grids_tmp%A1_gr,grids_tmp%rho_gr,grids_tmp%t_gr)
        else
            !Default guess. For consistency, send it to other images the same way as a guess
            !read from file. For small grids it takes a fraction of a second, for large grids
            !this will virtually never be used because we will have a value function from earlier
            !stage of CTS algorithm or a file.

            shared_data_V%V_tmp = 0.0_dp !(EFE) can find better constant guess
        end if
    else !(CTS_ind > 1)
        !load value function from previous stage of CTS algorithm using interpolation.
        !This is pretty much the same as when loading it from file, we just get it
        !from somewhere else.

        !We should already have allocated V_vect_tmp and grids_tmp from
        !previous stage of CTS algorithm and these
        !should contain values. If not, an error will result!
        if(.not. all([allocated(V_vect_tmp),allocated(grids_tmp%a1_gr),&
        allocated(grids_tmp%rho_gr),allocated(grids_tmp%t_gr)])) then
            call double_output(folder_name,'Error: when trying to read data from last stage of &
            CTS algorithm, something is not allocated.')
            error stop
        end if

       call V_from_V(reshape(V_vect_tmp,[N_a_tmp,N_rho_tmp,N_t_tmp,M_par]),&
        grids_tmp,grids,shared_data_V%V_tmp,folder_name)
        !We already checked that these variables are allocated
        deallocate(V_vect_tmp)
        deallocate(grids_tmp%a1_gr)
        deallocate(grids_tmp%rho_gr)
        deallocate(grids_tmp%t_gr)
    end if
!
end if !if(i_img == spec_img)
sync all !Make sure that the data is already available on special image before
!other images try to read it...


!V_old keeps, at every image, the copy of the entire value function from previous iteration.
allocate(V_old(N_a,N_rho,N_t,M_par))

!Retrieve the value function by all images and free up memory on special image
call get_V_from_spec(shared_data_V,V_old,spec_img)
sync all !have to sync before deallocating to make sure that the data is
!still there when other images try to read it!

!Now associate pointer with V_old (so that 'lower' subroutines can access it,
!the pointer is passed using common block). TBH could associate the pointer
!earlier and only once. Change at some point...
V_old_pntr => V_old


!if sim_only, then exit (this will exit the loop over CTS_ind...)
if(pars%sim_only) then
    exit
end if


!At this point, all images should have the value function and
!we can start VFI algorithm.

!Allocate memory to the stub of the value function in which every image will save
!results. (note that indices don't run from 1 but correspond
!to the indices in the whole unfolded value function).
allocate(V_stub(min_ind_img:max_ind_img))

!Allocate memory to policy functions
allocate(C_pol_stub(M_par,I_par,min_ind_img:max_ind_img),C_pol_stub_old(M_par,I_par,min_ind_img:max_ind_img))

!____Value function iteration__________________________


if(i_img == spec_img) then
    call double_output(folder_name,'Starting Value Function Iteration')
    VFI_stop_now = .false.
end if


!This is the variable to which all the images save their results.
!It should be allocated all the time on image spec_img, unallocated on all other images.
!(to save memory)
if(i_img == spec_img) then
    allocate(shared_data_V%V_stub_stacked(N_a*N_rho*N_t*M_par))
end if
sync all !Make sure that no image tries to save data there before it is allocated


do V_ind = 1,pars%VFI_max_iter
    !Cycle over gridpoints which this image works on
    !In the first n iterations (n should be small due to computational feasibility),
    !global optimization subroutines are possibly used - there are sometimes issues with
    !these in the case of multidimensional maximization (M>1) so need to be careful.

    if(V_ind<=pars%VFI_glob_opt_iter) then
        glob_opt = .true.
    else
        glob_opt = .false.
    end if
    loss_img_tot = 0.0_dp !initialize total loss at image

    !If we are skipping maximization step (Howard), set the value of skip_max accordingly,
    !and write a message if we are at special image
    if(skip_count_How > 0) then
        skip_max = .true. !change this to true later (stub so far)
        skip_count_How = skip_count_How -1
        if(i_img == spec_img) then
            call double_output(folder_name,'Skipping maximization step.')
            write(tmp_string1,'(I5)') skip_count_how
            call double_output(folder_name,trim(tmp_string1))
        end if
    else
        skip_max = .false.
    end if


    do gp_ind = min_ind_img,max_ind_img
        !The following subroutine solves the maximization problem at a single gridpoint

        call VFI_solve_gp(ind_unfold_all(:,gp_ind),grids,pars,N_a,N_rho,N_t,&
        c_guess_img(gp_ind),V_new_gp,skip_max,C_pol_gp,glob_opt,loss_gp,optim_wb_img(gp_ind))


        V_stub(gp_ind) = V_new_gp

        !Save the new policy function into the stub C_pol_stub (maybe call it C_pol_img)
        !and compute some statistics from it (like stopping criterion - the same as in case
        !of value function - and pass this across images like loss function)


        C_pol_stub(:,:,gp_ind) = C_pol_gp

        !update total loss at image
        loss_img_tot = loss_img_tot + loss_gp

    end do



    !After computing its part of value function, every image sends its peace of
    !it to the special image
    call send_V_to_spec(shared_data_V,V_stub,spec_img,min_ind_img,max_ind_img)
    sync all !Sync to make sure all pieces arrived

    !Compute difference in policy functions between iterations,
    C_pol_MARD_img = maxval(reshape(abs(C_pol_stub-C_pol_stub_old)/(abs(C_pol_stub_old) + 1.0_dp),&
    [M_par*I_par*(max_ind_img-min_ind_img+1)]))
    C_pol_AARD_img = sum(reshape(((abs(C_pol_stub-C_pol_stub_old)/(abs(C_pol_stub_old) + 1.0_dp))/&
    (M_par*I_par)),&
    [M_par*I_par*(max_ind_img-min_ind_img+1)]))

    sync all !Sync before recovering these numbers at the special image
    if(i_img == spec_img) then
        C_pol_MARD = 0.0_dp
        C_pol_AARD = 0.0_dp
        do img_ind = 1,n_img
            !Maximum absolute (value of) relative difference
            C_pol_MARD = max(C_pol_MARD,C_pol_MARD_img[img_ind])
            !Average absolute (value of) relative difference
            C_pol_AARD = C_pol_AARD + C_pol_AARD_img[img_ind]/(V_num_el)
        end do
    end if
    !The current policy function will be the 'old' one the next iteration.
    C_pol_stub_old = C_pol_stub

    !Compute average loss across all grid points.
    !loss_img_tot, loss_all_avg
    loss_all_avg = 0.0_dp
    do img_ind = 1,n_img
        loss_all_avg = loss_all_avg + loss_img_tot[img_ind]/(N_a*N_rho*N_t*M_par)
    end do


    !Now image spec_img has all the data. We can't read them directly (at other imagers) because of segmentation fault
    !errors (at least in the version of Intel compiler which I used - 15.0.something,
    !it could also be system-dependent), so
    !on image spec_img we first reshape the data and save the result in a temporary variable to which we allocate
    !memory. Then all other images retrieve data.
    if(i_img == spec_img) then
        if(.not. allocated(shared_data_V%V_tmp)) then
                 allocate(shared_data_V%V_tmp(N_a,N_rho,N_t,M_par))
        end if
        call V_unfold(shared_data_V)
    end if

    !Now that we have the new value function on image spec_img, compute some statistics
    !Maximum absolute relative difference b/w new and old value function
    if(i_img == spec_img) then
    V_MARD = maxval(reshape(abs((shared_data_V%V_tmp - V_old)/(abs(V_old)+1.0_dp)),[N_a*N_rho*N_t*M_par]))
    V_AARD = sum(reshape((abs(shared_data_V%V_tmp - V_old)/(abs(V_old)+1.0_dp))/&
        (N_a*N_rho*N_t*M_par),[N_a*N_rho*N_t*M_par]))
    end if

    !Report some stuff (later on write a proper subroutine for this).
    if(i_img == spec_img) then
        write(output_string,*) 'Iteration', V_ind
        call double_output(folder_name,trim(output_string))
        write(output_string,*) 'V_MaxARD = ',V_MARD
        call double_output(folder_name,trim(output_string))
        write(output_string,*) 'V_AvgARD = ',V_AARD
        call double_output(folder_name,trim(output_string))
        write(output_string,*) 'C_pol_MaxARD = ',C_pol_MARD,'C_pol_AvgARD = ',C_pol_AARD
        call double_output(folder_name,trim(output_string))
        !Also would be nice to compute average loss (for debugging purposes)
        if(.not. glob_opt) then !loss_gp in VFI_solve_gp is not computed correctly in this case
        write(output_string,*) 'loss_all_avg = ',loss_all_avg
        call double_output(folder_name,trim(output_string))
        end if
    end if

    !Stopping rule
    !(check only when we're not skipping maximization in Howard algorithm)
    if(i_img==spec_img .and. skip_max == .false.) then
        if(pars%stopping_rule_quantile > 0.0_dp .and. pars%stopping_rule_quantile < 1.0_dp) then
            !Quantile stopping rule
            !Compute the quantile (this can take a while so only do it if quantile stopping rule used)
            allocate(ARD_vector(N_a*N_rho*N_t*M_par))
            ARD_vector=reshape(abs((shared_data_V%V_tmp - V_old)/(abs(V_old)+1.0_dp)),[N_a*N_rho*N_t*M_par])
            ifail = 0
            quant_nag_input(1) = pars%stopping_rule_quantile
            call g01amf(N_a*N_rho*N_t*M_par,ARD_vector,1,quant_nag_input,V_quant_ARD,ifail)
            write(tmp_string1,'(F5.2)') pars%stopping_rule_quantile
            write(tmp_string2,*) V_quant_ARD(1)
            output_string = trim(tmp_string1)//' V_ARD quantile = '//trim(tmp_string2)
            call double_output(folder_name,trim(output_string))

            if(V_quant_ARD(1)<pars%stopping_rule_threshold) then
                VFI_stop_now = .true.
                output_string = 'Stopping rule ('//trim(tmp_string1)//' quant abs rel diff) satisfied, exiting VFI loop.'
                call double_output(folder_name,trim(output_string))
            end if

            deallocate(ARD_vector)
        else if(pars%stopping_rule_quantile == 1.0_dp) then
            !If V_MARD is less then stopping_rule threshold, stop now.
            if(V_MARD<pars%stopping_rule_threshold) then
                VFI_stop_now = .true.
                call double_output(folder_name,'Stopping rule (max abs rel diff) satisfied, exiting VFI loop.')
            end if
        else
            call double_output(folder_name,'Error: value of stopping_rule_quantile is wrong!')
            error stop
        end if
        !also check runtime limit
        if(time_real_acc() > runtime_limit_unix) then
            VFI_stop_now = .true.
            call double_output(folder_name,'VFI runtime limit reached, exiting VFI loop.')
        else if(time_real_acc() > runtime_limit_unix_CTS .and. CTS_ind < pars%CTS_split_times + 1) then
            !runtime_limit_unix_CTS <= runtime_limit_unix_CTS. This serves to prevent cases
            !where there is no convergence in first stages of CTS algorithm and we would get
            !stuck there for thousands of iterations. Before introducing Howard algorithm the
            !limit on number of iterations was rather low, now we can't use that one so we
            !use runtime limit and very high limit on number of iterations (because with Howard there
            !can be very many)
            VFI_stop_now = .true.
            call double_output(folder_name,'VFI (CTS) runtime limit reached, exiting VFI loop.')
        end if


    end if

    !If stopping criterion not satisfied, check whether we should skip some maximization steps (Howard)
    !Only check this if we're not in a stage where we are skipping maximization steps.
    if(i_img == spec_img .and. VFI_stop_now == .false. .and. use_How == .true. .and. skip_max == .false.) then
        !Depending on value of pars%crit_How, check the relevant criterion and decide
        !whether we're skipping maximization.

        select case(pars%crit_How)
            case(1)
                crit_value_How = V_MARD
            case(2)
                crit_value_How = V_AARD
            case(3)
                crit_value_How = c_pol_MARD
            case(4)
                crit_value_How = c_pol_AARD
            case default
                if(i_img == 1) then
                    write(*,*) 'Error: wrong value of parameter crit_value_How in the parameter file.'
                end if
                sync all
                error stop
        end select

        !If the critical value is lower than c_How, don't use Howard algorithm anymore. This is to
        !avoid overshooting.
        if(crit_value_How < pars%c_How) then
            use_How = .false.
            go to 666
        end if

        !Find the lowest i s.t. the criterion is lower than pars%a_How(i)
        do i_How = N_how,1,-1
            if(crit_value_How < pars%a_How(1,i_How)) then
                skip_count_How = pars%b_How(1,i_How)

                write(tmp_string1,*) '******* Howard acceleration algorithm criterion satisfied *******'
                write(tmp_string3,'(I5)') skip_count_HoW
                write(tmp_string2,*) 'Skipping ', trim(tmp_string3), ' maximization steps.'

                call double_output(folder_name,trim(tmp_string1))
                call double_output(folder_name,trim(tmp_string2))
                call double_output(folder_name,'*****************************************************')

                exit !no need to check other values as we assume that a(i) is increasing
                !and we care about the 'strictest' comparison only.
            end if
            !If this is never satisfied, the skip count will remain at zero (unless there is
            !a bug somewhere, this branch of program should only be reached with value 0)
        end do

666    end if
!pass the outcome to all images
sync all
use_How = use_How[spec_img]
skip_count_How = skip_count_How[spec_img]


    !Save value function (only if we're at the last stage of CTS, nobody cares about intermediate results
    !on coarse grid). Also don't save if save_iter_multiple is 0.
    if((mod(V_ind,pars%save_iter_multiple) == 0 .or. V_ind == pars%VFI_max_iter .or. VFI_stop_now)&
     .and. i_img==spec_img .and. CTS_ind == pars%CTS_split_times + 1 .and. pars%save_iter_multiple/=0) then
        call save_V(shared_data_V%V_stub_stacked,grids,N_a,N_rho,N_t,M_par,folder_name)
    end if

    if(i_img == spec_img) then
        !runtime report
        call runtime_report(initial_time,.true.,elapsed_time_string)
        call double_output(folder_name,trim(elapsed_time_string))
        call double_output(folder_name,'___________________________________________')
    end if

    sync all !So no image tries to read before 'unfolding' happens on special image

if(i_img==spec_img .and. pars%debug_mode == .true.) then
        call runtime_report(initial_time,.true.,elapsed_time_string)
        write(*,*) 'Runtime before pass = ',trim(elapsed_time_string)
end if

!debug
!write(*,*) 'debug: V_old(1,1,1,1) before get_V_from_spec',V_old(1,1,1,1)

!if(i_img == spec_img) write(*,*) 'debug: V_tmp(1,1,1,1) before get_V_from_spec',shared_data_V%V_tmp(1,1,1,1)
    call get_V_from_spec(shared_data_V,V_old,spec_img)

!write(*,*) 'debug: V_old(1,1,1,1) after get_V_from_spec',V_old(1,1,1,1)

    sync all !just in case we got to deallocation of shared_data_V%V_tmp
    !before all images finished getting the data - unlikely but not impossible


if(i_img == spec_img .and. pars%debug_mode == .true.) then
        call runtime_report(initial_time,.true.,elapsed_time_string)
        write(*,*) 'Runtime after pass = ',trim(elapsed_time_string)
end if


    !If VFI_stop_now == 1 at image spec_img, then stop. It is important that
    !there is a sync all statement b/w here and the point where the value is
    !possibly changed.
    if(VFI_stop_now[spec_img]) exit !
    !After this exit, all images should have access to the latest value function.


end do
!!_____________end of VFI_______________________________________

!Before the next iteration (or at the end if this is the end indeed), deallocate V_old (on image one, save
!it to a temporary variable if this is not the last iteration)

!deallocate these before going to the next state of CTS algorithm - the dimensions will be different in the next one.
deallocate(c_guess_img)
deallocate(optim_wb_img)
deallocate(V_stub)
!If this is the last iteration of CTS algorithm, keep V_old
if(CTS_ind < pars%CTS_split_times + 1) then
    deallocate(V_old)
end if
deallocate(ind_unfold_all)
deallocate(C_pol_stub,C_pol_stub_old)

!This thing could in theory be deallocated in every iteration over V_ind after
!passing data to free up memory - but in practice this leads to errors in passing
!data between images on different nodes when the passed arrays are very large
!(order of hundreds of MBs and more)!!! Some bug in CAF implementation...
if(i_img == spec_img .and. allocated(shared_data_V%V_tmp)) deallocate(shared_data_V%V_tmp)

!If this is not the last stage of CTS algorithm, save the value function in vector
!form and the old grids. Do this on image spec_img only.
if(i_img == spec_img .and. (CTS_ind < pars%CTS_split_times + 1)) then
    !save the old number of gridpoints (unnecessary here but used elsewhere so it's
    !easier and it doesn't matter
    N_a_tmp = N_a
    N_rho_tmp = N_rho
    N_t_tmp = N_t
    allocate(V_vect_tmp(N_a_tmp*N_rho_tmp*N_t_tmp*M_par))
    allocate(grids_tmp%A1_gr(N_a_tmp),grids_tmp%rho_gr(N_rho_tmp),&
        grids_tmp%t_gr(N_t_tmp))
    V_vect_tmp = shared_data_V%V_stub_stacked
    grids_tmp = grids
end if

!Once saved, deallocate the old grids. This sync all should actually not be necessary
!because here, image spec_img only reads its own copy of the data.
sync all
if(allocated(shared_data_V%V_stub_stacked)) deallocate(shared_data_V%V_stub_stacked)

!Keep the grids if this is the last stage of CTS
if(CTS_ind < pars%CTS_split_times + 1) then
    deallocate(grids%A1_gr,grids%rho_gr,grids%t_gr)
end if

end do !End of CTS algorithm (and thus VFI).

!sync after we got solution of the VFI problem
sync all

!______________________SIMULATION______________________________________________

!Here call subroutine sim_series, which will (for given value function
!and parameters, generate a simulated series of model variables).
if(i_img == spec_img) then
    write(*,*) 'Computing simulated series.'

    !Allocate memory to simulated series
    call sim_result_alloc(sim_result,pars,1)

    !and to the average series
    call sim_result_alloc(sim_avg,pars,1)

!Allocate memory to the type containing all of the simulations:
allocate(sim_res_all%SRA(pars%N_sim))
allocate(sim_res_all%tau_resid_last(pars%N_sim))
allocate(sim_res_all%tau_resid_MAD(pars%N_sim))
allocate(sim_res_all%tau_resid_var(pars%N_sim))
!Tracking quantiles of state variables (across simulations).
allocate(sim_res_all%rho_prime_05(pars%T_sim))
allocate(sim_res_all%rho_prime_95(pars%T_sim))
allocate(sim_res_all%a_prime_05(pars%T_sim))
allocate(sim_res_all%a_prime_95(pars%T_sim))
!and quantiles of tax policy
allocate(sim_res_all%tau_05(pars%T_sim))
allocate(sim_res_all%tau_95(pars%T_sim))
allocate(sim_res_all%trans_05(pars%T_sim))
allocate(sim_res_all%trans_95(pars%T_sim))


    !Solve first-period problem for every shock realization.
    do m_ind = 1,M_par
        call solve_FP_problem(V_old,pars%b_init,pars%theta_0,pars,FP_result(m_ind),grids%t_gr,m_ind)
    end do

do sim_ind = 1,pars%N_sim

    !Just for tracking progress in long simulations
    write(*,*) 'Computing simulation number ',sim_ind

    !Allocate memory to simulated series
    call sim_result_alloc(sim_res_all%SRA(sim_ind),pars,1)

    !Initialize RNG at some value - if we do it within subroutine that generates Markov chain,
    !there are segmentation faults - some issue with CAF which I couldn't figure out...
    !call random_seed(rng_seed) !use call random_seed without arguments for random initialization
    call random_seed


    !Use first-period solution as initial guess in the first maximization
    !(use the same value for all states if shocks are small), and
    !at this stage implement interior solution only. For corner solution we would
    !first need to look for feasible choices (otherwise we could get crazy results)
    !This would be done here - before calling subroutine sim_series - basically
    !just like LFFC but done at one point only so it should be quick.

    !subroutine sim_series should only generate the results and save them in a
    !defined type sim_result (like FP_result). And then another subroutine will
    !be implemented to save them (and plotting will then be done in Matlab).


    !In the 'first period', randomize over initial shock (unless initial shock given). Then pick the appropriate element
    !of FP_result.
    if(pars%s_init_index == -1) then
!                k,r1,r2
        if(M_par == 1) then
            s_init_ind_a = 1
        elseif(M_par == 2) then
            k = pars%P(M_par,1)/(1.0_dp + pars%P(M_par,1) - pars%P(1,1))
            !Just for debugging purposes compute the whole stationary probability vector (but k is all we need)
            !r1 = k
            !r2 = k * (1-pars%P(1,1))/pars%P(2,1)

            !k is the asymptotic (stationary) probability of state 1. So draw a number from U(0,1),
            !and if it's less than k, the initial state will be set to 1, otherwise it will be set to two.
            call random_number(randnum)
            if(randnum < k) then
                s_init_ind_a = 1
            else
                s_init_ind_a = 2
            end if
        else
            write(*,*) 'Error: need to generalize subroutine solve_FP_problem for M_par>2.'
            error stop
        end if
    else
        s_init_ind_a = pars%s_init_index
    end if

    !Use the first-period result correponding to the drawn shock.
    call sim_series(sim_result,FP_result(s_init_ind_a),pars,grids)

    !Now, if everything is OK, the whole sim_result should be possible to be saved at once.

    sim_res_all%SRA(sim_ind) = sim_result
end do

!If N_sim > 1 (compute statistics) of interest, average realization, residuals, ...
if(pars%N_sim > 1) then
    call sim_stat(sim_res_all,pars,sim_avg)
else
    call sim_stat(sim_res_all,pars,sim_avg)
    !write(*,*) 'N_sim = 1, therefore no statistics computed.'
end if

call save_stat(sim_res_all,folder_name,pars)

    !Save the sample average of simulated series. To save a particular
    !realization, set N_sim = 1 in parameter file.
    call save_series(sim_avg,folder_name,pars)

    !deallocate memory allocated to simulated series
    call sim_result_alloc(sim_result,pars,0)
    call sim_result_alloc(sim_avg,pars,0)
    deallocate(sim_res_all%SRA)

    deallocate(sim_res_all%tau_resid_last)
    deallocate(sim_res_all%tau_resid_MAD)
    deallocate(sim_res_all%tau_resid_var)

    deallocate(sim_res_all%rho_prime_05)
    deallocate(sim_res_all%rho_prime_95)
    deallocate(sim_res_all%a_prime_05)
    deallocate(sim_res_all%a_prime_95)

    deallocate(sim_res_all%tau_05)
    deallocate(sim_res_all%tau_95)
    deallocate(sim_res_all%trans_05)
    deallocate(sim_res_all%trans_95)
end if


!Display runtime and save it in log file
sync all
if (i_img == spec_img) then
    call runtime_report(initial_time,.true.,elapsed_time_string)
    call double_output(folder_name,trim(elapsed_time_string))
end if

end program taxhetgr
