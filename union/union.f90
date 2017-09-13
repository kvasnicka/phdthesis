program union
!This is the main file of the program taxhetgr (taxation under heterogeneous growth)


use mod_types !types for integer/real variables
use mod_IO !subroutines used for input/output
use mod_par !definition of type par and subroutine for loading this from parameter file
use mod_tools !various tools (tracking runtime, etc.)
use mod_parameters_hard !parameters which are hard-coded and not loaded from parameter file
use mod_work_dist !tools for distributing workload between images (processes)
use mod_union !types, subroutines, functions useful in this project only
use mod_utility

use nag_library

!The following module is for 1-6D bspline.
use bspline_sub_module


use ifport !so we can use system()

implicit none

!___________________________________________________________________
!Variable declaration: later on either put this into an include statement or a module
!- if it's to be used in other modules - probably rather use include statement.
!(because the set of variables used in various modules will not be the same)
character(256) :: parameter_file = ''

character(256) :: folder_name = ''
character(256) :: elapsed_time_string
type(par), codimension[*] :: pars !Type which contains all parameters, defined in module mod_par
!this can either be a coarray, or reading of data from parameter file could be done
type(par), target :: pars_img !Local copy of parameters, used because we need to associate its elements
!to a pointer and we can't give the TARGET attribute to a coarray.

real(dp) :: initial_time !Initial unix time including miliseconds.

!CAF control variables
integer :: i_img = 0
integer :: n_img = 0
integer :: V_num_el = 0 !number of elements of the value function.
integer :: min_ind_img = 0
integer :: max_ind_img = 0

!Value function
real(dp), allocatable, dimension(:,:,:,:,:), target :: V_old !Allocated all the time at all images.
real(dp), allocatable, dimension(:) :: V_stub !Piece of value function unfolded into row vector, allocated on all images

!also a copy for checking convergence if bsplines are used
real(dp), allocatable, dimension(:,:,:,:,:) :: V_old_copy


type(shared_data_V_type), codimension[*] :: shared_data_V
!allocated on special image (index spec_img) only, used for passing results
!(more precisely, this is allocated on all images, but its components are allocated to size 1
!on all images but one). Also note that reallocating components of this doesn't lead to implicit
!synchronization as reallocating whole coarrays does -> need to sync all manually!!!

type(shared_data_pol_type), codimension[*] :: pol !policy (both for sharing between images, and for keeping it)
!to replace x_pol_img eventually.

!Some VFI algorithm iteration-invariant data saved on every image:
integer, allocatable, dimension(:,:) :: ind_unfold_all !contains all of the 'unfolded' indices for every image,
!where unfolded index refers to the value function in the multi-dimensional as opposed the vector form.

logical, codimension[*] :: VFI_stop_now !This value will be determined on image spec_img only
!all images will stop if it is .true.

logical :: stop_now !this will be true on image 1 when file stop_now.txt is present in the results/folder_name folder

!CTS algorithm variables
integer :: CTS_ind = 0

real(dp), codimension[*] :: share_feas_img !weighted share of grid points at which feasible solution was found.
real(dp), codimension[*] :: avg_loss_img !average loss at image (loss from violating constraints)
!integer :: img_ind !used when we need to cycle over images

!variables for generating output (useful when we need to save formatted output).
character(256) :: output_string, tmp_string1,tmp_string2,tmp_string3

!Numbers of grid points (these are not the same as those in the parameter file because they can be
!changed by CTS algorithm)
integer :: N_a = 0 !number of grid points for asset holdings
integer :: N_rho = 0 !number of grid points for ratio of marginal utilities
integer :: N_aex = 0 !number of grid points for external debt
!(number of grid points per shock process is given in parameters_hard: M_par)

!Temporary variables: used for loading the value function from file and in between
!iterations of CTS algorithm.
integer :: N_a_tmp,N_rho_tmp,N_aex_tmp
real(dp), dimension(:), allocatable :: V_vect_tmp !temporary value function in vector form (used when loading
!value function from file and in CTS algorithm).

integer(dp) :: attmax_backup

real(dp) :: share_feas,avg_loss

integer :: V_ind !VFI algorithm index

type(grids_type), target :: grids, grids_tmp !grids and temporary grids (the latter used in CTS algorithm and
!when loading a value function from file.

type(ss_type), dimension(0:M_par) :: SS !steady-state. The first element (index 0) corresponds to steady state when
!government expenditure is kept forever at the expected (long-run SS probabilities) government expenditure.
!Element with index m contains the same, but the expenditure corresponds to shock realization with index m.
!This is useful for getting initial guess in first period, and for debugging.

!index for cycling over grid points
integer :: gp_ind

integer :: shep_ind

!index for cycling over images
integer :: img_ind

!Howard improvement algorithm
logical :: skip_max = .false. !If this is true, the maximization step will be skipped
!in the current iteration of VFi algorithm
real(dp) :: V_new_gp !new value of V at a grid point

integer :: rng_seed = 0 !some value for initializing RNG

!runtime_limit_unix: unix time at which VFI should stop
real(dp), codimension[*] :: runtime_limit_unix
real(dp), codimension[*] :: runtime_limit_unix_CTS !runtime limit for CTS algorithm
real(dp), codimension[*] :: runtime_limit_unix_How !runtime limit for Howard's algorithm.
real(dp)  :: time_unix
logical :: skip_check !if true, then stopping rule will not be checked.

!Pointer to value function. This will be used to make the old value function V_old
!accessible to the subroutine used in NAG maximization.
real(dp), dimension(:,:,:,:,:), pointer :: V_old_pntr
!Pointers for passing grids. Unfortunately, the type grids_type contains
!allocatable arrays and pointers to the whole type could not be part of a
!common block. Therefore I define pointers for each component of the grid.
!When we generalize to I>2, this will have to be changed everywhere (because
!the grids form will be different). Also create pointers for other things
!which need to be passed to some subroutines deep in the program.
real(dp), dimension(:), pointer :: a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr

real(dp), dimension(:), allocatable :: x_AAD_img !this is used only if rhoendadapt = .true. for dynamically changing
!he value of rhoend (if optim_subroutine == 1).
real(dp) :: x_AAD_placeholder = 1.0 !This is passed to subroutine VFI_solve_gp instead
!of an element of x_AAD_img if this is not allocated. Just to avoid some errors.

type(par), pointer :: pars_pntr !This will be viable only as long as pars contain
!no allocatable arrays. Right now it is easier to do this than reweite get_c_im
!and the function which it uses but when I generalize to I>2 and if I want the code
!to work for arbitrary M w/o recompilation, I will need to change this.
!___________________________________________________________________________

!Maximum (Average) Relative Difference in policy function between
!iterations. The average is also weighted by the number of grid points.
real(dp), codimension[*] :: x_MARD_img, x_AARD_img
real(dp), codimension[*] :: x_MARD, x_AARD
!Components for computing the above numbers (for each grid point)
real(dp) :: x_AARD_gp,x_MARD_gp
!variables for debugging
real(dp) :: debug_val,debug_val2

real(dp) :: loss_gp, loss_all_avg !loss at grid point, average loss at all grid points.
real(dp), codimension[*] :: loss_img_tot !loss average across gp's that an image works on

real(dp) :: V_MARD,V_AARD !maximum and average absolute relative difference.
real(dp), dimension(1) :: V_quant_ARD, quant_nag_input
real(dp), dimension(:), allocatable :: ARD_vector !vector of all absolute relative differences
integer :: ifail !for use in NAg
!_____________________________________________________________________
type(FP_result_type) :: FP_result !result of the first-period problem
type(FP_result_type), dimension(M_par) :: FP_result_allm !result of the first-period problem for
!all shocks realization. These are used in case of perm_eps = 0.0_dp (we can then save some
!time by computing the fist-period problem once only)
real(dp), dimension(ncvar_par) :: x_guess_t1 !initial guess of choice variables for period t=1.

type(sim_result_type) :: sim_result,sim_avg !simulated series and sample average
!_____________________________________________________________________

logical, codimension[*] :: use_How !is true if using Howard acceleration algorithm. Will be set to .false.
!once a criterion is less than c_How (see parameter file for definition and discussion)
integer, codimension[*] :: skip_count_How !counter of maximization steps to be skipped
real(dp) :: crit_value_How !value of the criterion
integer :: i_How

type(sim_result_all_type) :: sim_res_all !This contains all the simulated series (will be allocated on image 1 only),
!and some other results.
type(sim_result_all_type) :: sim_res_all_tmp !The same but temporary storage used for discarding outlying samples.
integer :: sim_ind = 1 !index for cycling over simulated samples.

real(dp), dimension(1,I_par) :: b_init_copy !copy of b_init
real(dp) :: b_ex_init_copy !copy of b_ex_init

integer, allocatable, dimension(:,:), target :: shep_inc !contains increments in the shepard's interpolation algorithm
!(see subroutine V_interp_shep in mod.union.f90 for how it's used). Essentially we precompute a small matrix
!so that we can do 1 do cycle deep in the program instead of 4 cycles.
!Also create a pointer so we can put shep_inc into common blocks
integer, pointer, dimension(:,:) :: shep_inc_pntr


integer :: disc_num !number of discarded samples.
real(dp), dimension(:), allocatable :: disc_crit !array containing value of discarding criterion for all samples
integer, dimension(:), allocatable :: disc_indices !indices of discarded samples
integer :: ind_free !contains the next available free index in the temporary storage

logical :: file_exists

integer :: m_ind
integer :: s_init_ind_a = 1 !the actual initial period index (in case of randomization)
        !The following variables are use for randomly drawing initial shock realization
        !if s_init_index = -1 in parameter file. We want to draw the shock with frequency associated
        !with stationary distribution.
        real(dp) :: k,r1,r2,randnum

integer :: maxatt_iter !the maximum number of attempts in maximization at a given iteration in VFI algorithm.
real(dp) :: attmax_rndnum
logical :: attmax_special
logical :: in_the_box

integer, allocatable :: s_ind_loaded(:,:)


real(dp) :: V_avg !(used only if pars%debug_mode)

integer :: N_sim_done !the number of simulations which were actually done before runtime expired

logical :: DA_this_iter !a logical variable which will be true if deallocations are supposed to be made
!in an iteration of CTS. This is almost always true and it is false only if NDA = .true. in parameter file,
!and we are in the last iteration of CTS algorithm. The role of this is to prevent some weird bugs on HPC

!variables for PSO in VFI
logical :: PSO_this_iter = .false.
real(dp) :: PSO_rndnum = 0.0_dp
real(dp) :: VFI_rndnum = 0.0_dp

!debug variable
integer, dimension(:), allocatable :: CPO_optim_count

!============variables for bspline===============
integer :: bspline_iflag !gives some output info (if it is 0 on output no errors).
!k's are the order of splines (2 is linear, 4 is cubic).
!test is with 2 to begin with.

!t's are the arrays of knots
real(dp), dimension(:,:), allocatable, target :: t_a1,t_a2,t_rho,t_aex
!and pointers so we can put these vairbales in common block
real(dp), dimension(:,:), pointer :: t_a1_pntr,t_a2_pntr,t_rho_pntr,t_aex_pntr
integer, dimension(M_par) :: inbvx,inbvy,inbvz,inbvq,iloy,iloz,iloq
!initialization parameter which must be set to 1 the first time this routine is called
!(every time that I recompute the bspline coefficients).
integer :: bspline_k_a1_backup,bspline_k_a2_backup,bspline_k_rho_backup,bspline_k_aex_backup

!===========end of variables for testing bspline===============
integer :: seed_n !size od the seed
integer, allocatable :: seed_array(:)


!This is a coarray which contains all the data which are exchanged between images
!during sim_series (if I parallelize solve_FP_problem eventually, then I will need
!a different type with different dimensions of choice variables).
type(sim_shared_type), codimension[*] :: sim_shared


common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
!to some subroutines.
common /pars_pointer/ pars_pntr

common /grids/ a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr

common /shep_pointer/ shep_inc_pntr !pointer used to pass a pre-computed matrix of increments in Shepard's interp.

!pointer used to pass information used in bspline evaluation (pass this directly to the interpolating subroutine)
!and not to CPC_rets. I want these subroutines to have essentially the same interface as CPC_ret.
common /bspline_pointer1/ inbvx, inbvy, inbvz, inbvq, iloy, iloz, iloq,&
t_a1_pntr,t_a2_pntr,t_rho_pntr,t_aex_pntr
!___________________________________________________________________________

!This follows (roughly) the guide for CAF debugging for Windows but was adapdted for ubuntu.
!it needs to be commented out in the version for HPC (or the Makefile needs to be updated to include this).
!Source:
!https://software.intel.com/en-us/articles/how-to-debug-fortran-coarray-applications-on-windows

!call ECD !(enable coarray debugging)

!Include taxhetgr_initialize.f90, which contains tasks such as loading
!parameters from file, etc. Some variables are also defined there and pointers
!associated.
include 'include/taxhetgr_initialize.f90'
!______________________________________________________________________________


!If pars%shocks_gen == .true., then the program only generates the shocks, then exits.
if(pars%shocks_gen_only) then
    if(pars%rng_seed < 0) then
        call random_seed !if seed negative,just initialize RNG randomly.
    else
        call random_seed(pars%rng_seed) !initialize rng seed for repeatability
    end if
    if(this_image() == 1) then
        call double_output(folder_name,'Shocks_gen is .true. in parameter file.')
        call double_output(folder_name,'The program will generate shocks and exit.')
        call double_output(folder_name,'These are saved in the same results folder as standard output.')
        call shocks_gen_and_save(folder_name,pars)
        goto 1988 !go to the end of the program.
    else
        goto 1988
    end if
end if




!To avoid confusion with deallocation, allocate the following small array used in bspline
!(if bspline is used) only once, and then refer to parts of it in earlier stages of CTS algorithm.
if(pars%VFI_interpolation_mode == 3) then
    allocate(t_a1(pars%N_a + pars%bspline_k_a1,pars%M_grid),t_a2(pars%N_a + pars%bspline_k_a2,pars%M_grid),&
    t_rho(pars%N_rho + pars%bspline_k_rho,pars%M_grid),t_aex(pars%N_aex + pars%bspline_k_aex,pars%M_grid))
    !and initialize these
    t_a1 = 0.0_dp;t_a2 = 0.0_dp;t_rho = 0.0_dp;t_aex = 0.0_dp;
end if

!initialize backup (it is done like this because this is passed quite deep into the program
!so it's easier than explicitly passing it to all the subroutines below
bspline_k_a1_backup = pars%bspline_k_a1
bspline_k_a2_backup = pars%bspline_k_a2
bspline_k_rho_backup = pars%bspline_k_rho
bspline_k_aex_backup = pars%bspline_k_aex

do CTS_ind = 1,pars%CTS_split_times + 1
!This part - possibly move to a file which will the be included (to make the code less messy).

!Figure out whether we are doing deallocation in this iterations of CTS. The only case in which
!I don't want to deallocate is when we are in the last stage of CTS and pars%NDA is true (not deallocate).
!Then take a .not. and if the resulting variable is true, I want to deallocate.
!If sim_interp then also do not deallocate because I will need the policy function in simulation.
DA_this_iter = .not. ((pars%NDA .or. pars%sim_interp) .and. CTS_ind == pars%CTS_split_times + 1)

!Howard acceleration algorithm - at every stage of CTS initialize the variable telling the
!program whether to use acceleration. It might also be a good idea to do this at the last
!stage of CTS only - this condition can be implemented here
if(pars%crit_How /= 0) then
    use_How = .true.
else
    use_How = .false.
end if
skip_count_How = 0 !Counter of maximization steps skipped - will be set to some positive value
!whenever a criterion is satisfied and every time a maximization step is skipped this counter
!will decrease by one.
skip_max = .false.

!Compute number of grid points at this iteration of CTS algorithm (at every grid point). If CTS
!algorithm is not used, then the values in parameter file will be returned by the subroutine.
call CTS_num_points(CTS_ind,pars,N_a,N_rho,N_aex)

!If we want to get symmetric grid for rho (rho_gr_symm = .true. in parameter file), make sure
!that the number of grid points N_rho is odd. If it is not, then set it to odd and display a message.
if(pars%rho_gr_symm .and. (N_rho/2)*2 == N_rho) then
    if(i_img == 1) call double_output(folder_name,'Using symmetric grid for rho. Setting N_rho odd.')
    N_rho = N_rho + 1
end if

V_num_el = N_a*N_a*N_rho*N_aex*pars%M_grid !number of elements of value function (used often)
!Note that in case of IID shocks, M_grid is 1 (not M_par) - because the value is independent of shock
!realization.

!If the number of gridpoints is less than a given constant, use only linear splines
!in this stage of CTS (for that variable).
!Actually, for some reason this seems to be causing problems sometimes. It's better not to use it
!(and set the minimum to pars%bspline_N_nonlin to 1 (and perhapsvoid using CTS altogether).
if(pars%VFI_interpolation_mode == 3) then
    if(N_a < pars%bspline_N_nonlin) then
        pars%bspline_k_a1 = min(2,pars%bspline_k_a1)
        pars%bspline_k_a2 = min(2,pars%bspline_k_a2)
    else !otherwise restore the order from backup
        pars%bspline_k_a1 = bspline_k_a1_backup
        pars%bspline_k_a2 = bspline_k_a2_backup
    end if
    if(N_rho < pars%bspline_N_nonlin) then
        pars%bspline_k_rho = min(2,pars%bspline_k_rho)
    else !otherwise restore the order from backup
        pars%bspline_k_rho = bspline_k_rho_backup
    end if
    if(N_aex < pars%bspline_N_nonlin) then
        pars%bspline_k_aex = min(2,pars%bspline_k_aex)
    else !otherwise restore the order from backup
        pars%bspline_k_aex = bspline_k_aex_backup
    end if
end if


if(i_img == spec_img .and. .not.(pars%sim_only)) then !Write this even if CTs not used (then CTS_ind = 1 out of 1)
    write(tmp_string1,'(I3)') CTS_ind
    write(tmp_string2,'(I3)') pars%CTS_split_times + 1
    output_string = 'CTS algorithm iteration '//trim(tmp_string1) //' out of '//trim(tmp_string2)
    call double_output(folder_name,trim(output_string))

    write(tmp_string1,'(I5)') N_a
    write(tmp_string2,'(I5)') N_rho
    write(tmp_string3,'(I5)') N_aex

    output_string = 'Number of gridpoints: N_a = '//trim(tmp_string1)//', N_rho = '//trim(tmp_string2)// &
    ', N_aex = '//trim(tmp_string3)
    call double_output(folder_name,trim(output_string))

    if(pars%iid_shocks) then
        call double_output(folder_name,'Shocks are i.i.d., setting M = 1 in grid construction.')
    else
        write(tmp_string1,'(I5)') M_par
        output_string = 'Shocks are not i.i.d., setting M = '// trim(adjustl(tmp_string1)) //' in grid construction'
        call double_output(folder_name,trim(output_string))
    end if

    write(tmp_string1,'(I20)') V_num_el
    output_string = 'Total number of gridpoints = '//trim(adjustl(tmp_string1))
    call double_output(folder_name,trim(output_string))

    write(tmp_string1,'(I6)') n_img
    write(tmp_string2,'(F15.1)') real(V_num_el)/real(n_img)
    output_string = 'Working on '//trim(adjustl(tmp_string1))//' images ('//trim(adjustl(tmp_string2))//' gridpoints per image)'
    call double_output(folder_name,trim(output_string))
    call double_output(folder_name,'____________________________________________________')
end if

!Get the range of indices (in the vector representation of value function)
!on which every image will work.
call get_img_ind(i_img,n_img,V_num_el,min_ind_img,max_ind_img)

!allocate memory to the unfolded indices. (indexed the same way as V_stub, so that
!element with index i contains the unfolded index for i-th term in the vector representation
!of the value function.
allocate(ind_unfold_all(5,min_ind_img:max_ind_img)) !5 state variables -> This already assumes I_par=2

!Get the unfolded indices (indices in the multidimensional value function corresponding to the
!indices in vector form value function)
call get_ind_unfold(min_ind_img,max_ind_img,[N_a,N_a,N_rho,N_aex,pars%M_grid],ind_unfold_all)

!Generate grids:
allocate(grids%A1_gr(N_a),grids%A2_gr(N_a),grids%rho_gr(N_rho),grids%aex_gr(N_aex))

call grids_gen(grids%A1_gr,grids%A2_gr,grids%rho_gr,grids%aex_gr,N_a,N_rho,N_aex,pars)

!Associate pointers with elements of the grids
a1_gr_pntr => grids%A1_gr
a2_gr_pntr => grids%A2_gr
rho_gr_pntr => grids%rho_gr
aex_gr_pntr => grids%aex_gr

!And the same for pointers used in passing information to V_interp_spline
t_a1_pntr => t_a1
t_a2_pntr => t_a2
t_rho_pntr => t_rho
t_aex_pntr => t_aex

!Get deterministic steady-state associated with initial asset holdings.

!Right now the steady-state returns variables associated with a steady state with no debt in autarky
!This should be a pretty ok choice because this is just supposed to provide a feasible initial guess,
!It can be far from optimum. And we are permutating the guess anyway. It's hard to come up
!with something better.
!This is needed to get initial guess both in LFFC, and in solving the first period problem.
!therefore it needs to be before checking pars%sim_only which would then skip it.
call get_aut_ss(pars,SS(0),-1)
do m_ind = 1,M_par
    call get_aut_ss(pars,SS(m_ind),m_ind)
end do

!if we only want to simulate, skip LFFC (go to 999
!which is loading of the value function). But before that load the policy function on image 1
!if sim_interp = .true.
if(pars%sim_only) then


if(pars%sim_interp) then

if(n_img > 1 .and. this_image()==1) then
    write(*,*) 'When sim_interp = .true., we might run out of memory.'
    write(*,*) '(consider using command export FOR_COARRAY_NUM_IMAGES=1)'
end if

    inquire(file=('results/'//trim(pars%input_folder)//'/x_pol_all.out'), exist=file_exists)
    if (.not. file_exists) then
        if(i_img == spec_img) then
            call double_output(folder_name,'sim_interp true but file x_pol_all.out not found')
            call double_output(folder_name,'Setting sim_interp = .false.')
        end if
        pars%sim_interp = .false.
        goto 999 !go to loading the value function
    end if

!On image 1, allocate memory to policy function and load it from a file.
    if(i_img == spec_img .or. pars%NDA) then
        allocate(pol%x_pol_all(ncvar_par,V_num_el))
    else !this prevents crashes even though it should not be necessary (specific to compiler)
        allocate(pol%x_pol_all(1,1))
    end if

    sync all


    !image spec_img loads the data, the other images then recover the part relevant to gridpoints
    !on those images only (this is important to save memory because the policy function is even a larger
    !object than the value function so we often can't afford all images to have the copy of the whole
    !policy function (and they don't need it for anything).
    if(i_img == spec_img) then
        call load_N(N_a_tmp,N_rho_tmp,N_aex_tmp,pars%input_folder,folder_name,pars%M_grid)
        !if any of the dimensions are wrong, ignore the policy function and write a warning
        if(N_a_tmp /= N_a .or. N_rho_tmp /= N_rho .or. N_aex_tmp /= N_aex) then
            call double_output(folder_name,'Warning: initial guess of the wrong dimension!')
            call double_output(folder_name,'Only the initial guess of value function will be loaded.')
            pars%sim_interp = .false.
        end if
        call load_policy(pol%x_pol_all,pars%input_folder,folder_name)
    end if

    sync all
    pars = pars[spec_img] !just in case we set pars%sim_interp = .false. on image spec_img.

    !Now the policy is loaded on image 1 only, it is not sent to other images. But that is fine
    !since simulation is done on one image only, so the other images do not need the policy function.

end if !pars%sim_interp

goto 999 !go to loading value function.

end if

!__________________INITIAL GUESS of optimal consumption___________________________
!Before starting value function iteration, we need to get an initial guess at every grid point.
!To begin with, we look for any feasible choice. In VFI algorithm, the initial guess
!will be updated at every iteration so that we hopefully start close to the optimum and
!save time. In later stages of VFI algorithm, Howard improvement algorithm will be used
!which essentially means that a consumption choice obtained as optimum at an iteration
!of VFI algorithm will be used for several iterations (skipping the maximization step).

!Allocate memory to policy function. This is in a type pol (coarray), which contains an allocatable
!array. This is for sharing policy function between images (useful mainly in between CTS algorithm
!stages, and when we save/load policy function)
allocate(pol%x_pol_img(ncvar_par,min_ind_img:max_ind_img))

!If pars%load_initial_guess == .true., then we are going to load initial guess of value
!function from file. If the policy function x_vect.out is also present in the folder,
!we will use that as initial guess. The policy function is assumed to have
!been obtained on the same grid (the main application is continuing computation which
!was not finished for some reason, such as running out of time on HPC).
if(pars%load_initial_guess) then
    !pars%M_grid replacing M_par here due to IID shocks possibility

    !Check if the initial guess of consumption is present. If not, proceed directly to
    !end of this branch of the program an start with looking for initial guess by a grid search.
    inquire(file=('results/'//trim(pars%input_folder)//'/x_pol_all.out'), exist=file_exists)
    if (.not. file_exists) then
        if(i_img == spec_img) then
            call double_output(folder_name,'Load_initial_guess true but file x_pol_all.out not found')
            call double_output(folder_name,'Only the initial guess of value function will be loaded.')
        end if
        go to 234
    end if

    !allocate memory for loading the data
    !DEBUG: if pars%NDA, then we are not deallocating memory (further below). What's more
    !we wastefully allocate memory for the whole policy function on ALL images. I have a feeling that
    !This asymmetric allocation in an element of a coarray of types is what caused problems on HPC!
    !(on multiple nodes, not on 1 - i.e., when coarray compiler option is distributed).
    if(i_img == spec_img .or. pars%NDA) then
        allocate(pol%x_pol_all(ncvar_par,V_num_el))
    else !this prevents crashes even though it should not be necessary (specific to compiler)
        allocate(pol%x_pol_all(1,1))
    end if

    sync all

    !image spec_img loads the data, the other images then recover the part relevant to gridpoints
    !on those images only (this is important to save memory because the policy function is even a larger
    !object than the value function so we often can't afford all images to have the copy of the whole
    !policy function (and they don't need it for anything).
    if(i_img == spec_img) then
        call load_N(N_a_tmp,N_rho_tmp,N_aex_tmp,pars%input_folder,folder_name,pars%M_grid)
        !if any of the dimensions are wrong, ignore the policy function and write a warning
        if(N_a_tmp /= N_a .or. N_rho_tmp /= N_rho .or. N_aex_tmp /= N_aex) then
            call double_output(folder_name,'Warning: initial guess of the wrong dimension!')
            call double_output(folder_name,'Only the initial guess of value function will be loaded.')
            go to 234
        end if

        call load_policy(pol%x_pol_all,pars%input_folder,folder_name)

    end if
    sync all !all images must wait till image spec_img loads the data.

    !Copy the data from image 1 to all images
    pol%x_pol_img(1:ncvar_par,min_ind_img:max_ind_img) = pol[spec_img]%x_pol_all(1:ncvar_par,min_ind_img:max_ind_img)

    !free up memory (have to sync before that to make sure data was sent to all the images)
    sync all


    if(DA_this_iter) deallocate(pol%x_pol_all)
    sync all

    goto 999 !skipping the LFFC algorithm since we already have a guess of policy function loaded from a file.
234 end if

if(i_img == spec_img) then
    call double_output(folder_name,'Looking for feasible choices...')
end if

!If CTS_ind > 1, we need to have the old grids in grids_tmp. Check that these are allocated
if(.not. (allocated(grids_tmp%a1_gr) .and. allocated(grids_tmp%a2_gr)&
 .and. allocated(grids_tmp%rho_gr) .and. allocated(grids_tmp%aex_gr)) .and. CTS_ind>1) then
    write(*,*) 'Error: grids_tmp not allocated when CTS_ind>1 in taxhetgr.f90!'
    sync all
    error stop
end if


call get_guess_img(pol%x_pol_img,SS,min_ind_img,max_ind_img,ind_unfold_all,pars,grids,&
share_feas_img,avg_loss_img,CTS_ind,pol%x_pol_all_unf,grids_tmp)

!If we have the grids deallocate memory (only on images except for spec_img, where they
!will be needed (possibly) for value function interpolation.
if(i_img /= spec_img) then
    if(allocated(grids_tmp%A1_gr) .and. DA_this_iter) deallocate(grids_tmp%A1_gr)
    if(allocated(grids_tmp%A2_gr) .and. DA_this_iter) deallocate(grids_tmp%A2_gr)
    if(allocated(grids_tmp%rho_gr) .and. DA_this_iter) deallocate(grids_tmp%rho_gr)
    if(allocated(grids_tmp%aex_gr) .and. DA_this_iter) deallocate(grids_tmp%aex_gr)
end if

!if pol%x_pol_all was allocated from previous stage of CTS alorithm, get rid of it
if(allocated(pol%x_pol_all_unf)) then
    if(DA_this_iter) deallocate(pol%x_pol_all_unf)
end if

!Sync after every image gets an initial guess of consumption. This is mainly so that we
!can compute the overall feasibility share and display a message about it (every image
!only needs its own copy of c_guess_img, no need to access initial guesses for gridpoints
!which are processed by other images).
sync all

!Display a message about results of LFFC
call LFFC_msg(share_feas_img,avg_loss_img,min_ind_img,max_ind_img,spec_img,V_num_el,folder_name,initial_time)


!The following section serves to get an initial guess of value function to all images. If CTS
!algorithm is used, and we are at the first iteration, we either load the guess from file
!by one image and send it to all images, or we do the same thing with a constant guess.
!If we are not at the first stage of CTs algorithm, one image computes a value
!function using interpolation from the value function computed on the coarse grid. This
!function is then sent to all images.

!We prepare the data to be sent on one image, then we synchronize and send it to all images.
999 if(i_img==spec_img) then

    !Allocate memory to object where we will save the new value function on one image only
    allocate(shared_data_V%V_tmp(N_a,N_a,N_rho,N_aex,pars%M_grid))

    if(CTS_ind == 1) then
        !Either load an initial guess from input_folder (given in parameter file) or use an initial guess
        !given here (or in an include file if it's something more complicated)
        if(pars%load_initial_guess==.true. .and. (.not. pars%gen_transfers_ss)) then
            !load the guess from old file on special image
            !!!Remember to use pars%M_grid instead of M_par here when generalizing due to IID shocks...

            !Get number of grid points from file
            call load_N(N_a_tmp,N_rho_tmp,N_aex_tmp,pars%input_folder,folder_name,pars%M_grid)


            !allocate memory to the temporary value function and grids and load these from files
            allocate(V_vect_tmp(N_a_tmp*N_a_tmp*N_rho_tmp*N_aex_tmp*pars%M_grid))
            allocate(grids_tmp%A1_gr(N_a_tmp),grids_tmp%A2_gr(N_a_tmp),grids_tmp%rho_gr(N_rho_tmp),&
                grids_tmp%aex_gr(N_aex_tmp))

            call load_V(V_vect_tmp,grids_tmp%a1_gr,grids_tmp%a2_gr,grids_tmp%rho_gr, &
                grids_tmp%aex_gr,pars%input_folder,folder_name)


            !(EFE) could paralelize the interpolation but there is not much to gain,
            !the code spends very little time here compared to other places.

            !reshaping the vector form value function so that is can be an input for V_from_V
            !This depends on the order of gridpoints in the vector form being correct.

            call V_from_V(reshape(V_vect_tmp,[N_a_tmp,N_a_tmp,N_rho_tmp,N_aex_tmp,pars%M_grid]),&
                grids_tmp,grids,shared_data_V%V_tmp,folder_name)

            !Now the interpolated value function is saved in shared_data_V%V_tmp, which will
            !later be retrieved by all images.

            !We already have the new value function, do not need to keep the old value
            !function or the old grid
            if(DA_this_iter) deallocate(V_vect_tmp)
            if(DA_this_iter) deallocate(grids_tmp%A1_gr,grids_tmp%A2_gr,grids_tmp%rho_gr,grids_tmp%aex_gr)
        else
            !Default guess. For consistency, send it to other images the same way as a guess
            !read from file. For small grids it takes a fraction of a second, for large grids
            !this will virtually never be used because we will have a value function from earlier
            !stage of CTS algorithm or a file.

            shared_data_V%V_tmp = 0.0_dp !(EFE) can find better constant guess (maybe just steady-state utility
            if(this_image() == 1) then
                write(*,*) 'Using naive guess of value function V = 0.'
            end if

            !discounted).
        end if
    else !(CTS_ind > 1)
        !load value function from previous stage of CTS algorithm using interpolation.
        !This is pretty much the same as when loading it from file, we just get it
        !from somewhere else.

        !We should already have allocated V_vect_tmp and grids_tmp from
        !previous stage of CTS algorithm and these
        !should contain values. If not, an error will result!
        if(.not. all([allocated(V_vect_tmp),allocated(grids_tmp%a1_gr),allocated(grids_tmp%a2_gr),&
        allocated(grids_tmp%rho_gr),allocated(grids_tmp%aex_gr)])) then
            call double_output(folder_name,'Error: when trying to read data from last stage of &
            CTS algorithm, something is not allocated.')
            error stop
        end if


       call V_from_V(reshape(V_vect_tmp,[N_a_tmp,N_a_tmp,N_rho_tmp,N_aex_tmp,pars%M_grid]),&
        grids_tmp,grids,shared_data_V%V_tmp,folder_name)

        !We already checked that these variables are allocated
        if(DA_this_iter) deallocate(V_vect_tmp)
        if(DA_this_iter) deallocate(grids_tmp%a1_gr)
        if(DA_this_iter) deallocate(grids_tmp%a2_gr)
        if(DA_this_iter) deallocate(grids_tmp%rho_gr)
        if(DA_this_iter) deallocate(grids_tmp%aex_gr)
    end if

end if !if(i_img == spec_img)

sync all !Make sure that the data is already available on special image before
!other images try to read it...

!V_old keeps, at every image, the copy of the entire value function from previous iteration.
allocate(V_old(N_a,N_a,N_rho,N_aex,pars%M_grid))

!If bspline_check_V_convergence, then also allocate a space for copy of V_old on images 1 for chcking convergence
if(pars%bspline_check_V_convergence .and. this_image()==1 .and.  pars%VFI_interpolation_mode == 3) then
    if(allocated(V_old_copy)) deallocate(V_old_copy)
    allocate(V_old_copy(N_a,N_a,N_rho,N_aex,pars%M_grid))
end if


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
    !(IMP) should deallocate memory allocated to some variables which will not be deallocated
    exit
end if

!Don't deallocate. Technically it would be nice to save memory but
!it could lead to problems in passing V between images for high
!number of gridpoints... (it lead to problems in the past despite using sync
!all, etc. A compiler-specific bug).
!if(allocated(shared_data_V%V_tmp)) deallocate(shared_data_V%V_tmp)

!At this point, all images should have the value function and
!we can start VFI algorithm.

!Allocate memory to the stub of the value function in which every image will save
!results. (note that indices don't run from 1 but correspond
!to the indices in the whole unfolded value function).
allocate(V_stub(min_ind_img:max_ind_img))

!The following array is used only for adapting rhoend.
if(pars%rhoendadapt) then
    allocate(x_AAD_img(min_ind_img:max_ind_img))
    !IF we are in the last stage of CTS set it to E-05 (so that rhoend will never
    !be less than E-06 in the first iteration of CTS. This might be the case if we decide
    !at some point to increase rhoend quite a bit (perhaps to E-03) to save time in regions
    !of state space where it does not matter.
    if(CTS_ind == pars%CTS_split_times + 1) then
        x_AAD_img = 0.00001_dp
    else
        x_AAD_img = 1.0_dp !initialize (this has to be as great
        !or greater than pars%rhoend, otherwise this will cause issues in the
        !first iteration. Also this can not be negative.
    end if
end if


!____Value function iteration__________________________
if(i_img == spec_img) then
    call double_output(folder_name,'Starting Value Function Iteration')
    call double_output(folder_name,' ')
    VFI_stop_now = .false.
end if


!This is the variable to which all the images save their results.
!It should be allocated all the time on image spec_img, unallocated on all other images.
!(to save memory)
if(i_img == spec_img) then
    allocate(shared_data_V%V_stub_stacked(N_a*N_a*N_rho*N_aex*pars%M_grid))
end if
sync all !Make sure that no image tries to save data there before it is allocated


do V_ind = 1,pars%VFI_max_iter
    !Cycle over gridpoints which this image works on
    loss_img_tot = 0.0_dp !initialize total loss at image

    if(this_image() == 1) then
        write(output_string,*) 'Iteration', V_ind
        call double_output(folder_name,trim(output_string))
    end if

    !If VFI_interpolation_mode is 3, then compute the coefficients for bspline.
    !This overwrite the value of V_old which will now contain the coeffients rather than the values at grid points.
    !In the future (if this is prohibitively long) I can try to parallelize this. At the very least, since
    !the value function is treated separately for each shock realization, I can do this on M_par images, with
    !1 image doing the interpolation. But with the number of grid points which I used this step does not
    !take very long and it is so much easier to write it like this...
    if(pars%VFI_interpolation_mode == 3) then

    !create copy of V_old on image 1 for convergence purposes
    if(pars%bspline_check_V_convergence .and. this_image() == 1) then
        V_old_copy = V_old
    end if


!initialize the following to 1 (and do not change them again).
inbvx = 1
inbvy = 1
inbvz = 1
inbvq = 1
iloy = 1
iloz = 1
iloq = 1

do m_ind = 1,pars%M_grid
    call db4ink(grids%a1_gr,N_a,grids%a2_gr,N_a,grids%rho_gr,N_rho,grids%aex_gr,N_aex,&
    V_old(:,:,:,:,m_ind),pars%bspline_k_a1,pars%bspline_k_a2,pars%bspline_k_rho,pars%bspline_k_aex,0,&
    t_a1(1:(N_a + pars%bspline_k_a1),m_ind),t_a2(1:(N_a + pars%bspline_k_a2),m_ind),&
    t_rho(1:(N_rho + pars%bspline_k_rho),m_ind),t_aex(1:(N_aex + pars%bspline_k_aex),m_ind),&
    V_old(:,:,:,:,m_ind),bspline_iflag)
end do


end if

!____end of computing bspline coefficients....


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

    !Initialize statistics
    x_MARD_img = 0.0_dp ; x_AARD_img = 0.0_dp

    !Get the maximum number of attempts (this is always one unless we are in the first
    !iteration of VFI and CTS). Also, the maximum number of attempts will not be more than
    !1 if we are loading initial guess (that is used when we continue with a previously
    !obtained value and policy function - the number of GPs is usually great in this case and we cannot afford to do this.

    if((CTS_ind == 1 .and. V_ind == 1 .and. pars%load_initial_guess == .false.) .or. &
    (V_ind == 1 .and. pars%always_perm_first_iter) ) then
        maxatt_iter = max(pars%CTS_first_maxatt,1)
        !attmax_special tells the program that this iteration of VFI is treated separately,
        !so that the value of maxatt_iter should not be changed (no matter what the values
        !of VFI_attmax, and VFI_attmax_rnd.
        attmax_special = .true.
    else !In all other cases stick to the default value
        maxatt_iter = pars%VFI_attmax !this value will actually not be used but it's better to set it to avoid bugs
        attmax_special = .false.
    end if
    !cycle over all grid points handled by an image
    do gp_ind = min_ind_img,max_ind_img

    !If this iteration is not 'special case', determine the value of attmax (number of atempts
    !in maximization at this gridpoint). Otherwise the value should already be correct in maxatt_iter.
    if(.not. attmax_special) then
    select case(pars%VFI_attmax_rnd)
        case(0) !No randomization is performed
            maxatt_iter = pars%VFI_attmax
        case(1) !Randomization with the smae chance to perturbate initial guess at all GPs
            call random_number(attmax_rndnum)
            !The probability used in this case is VFI_attmax_p1.
            if(attmax_rndnum <= pars%VFI_attmax_p1) then
                maxatt_iter = pars%VFI_attmax
            else !the default case (1 attempt)
                maxatt_iter = 1
            end if
        case(2)
            maxatt_iter = pars%VFI_attmax
            !No randomization in this case (and it is handled a bit more further below).
        case default
            !The parameter value will be the same on all images so we can afford the sync all
            if(this_image() == 1) then
                write(*,*) 'Error: Wrong value of parameter VFI_attmax_rnd'
                error stop
            end if
            sync all
    end select
    end if



    !This is rough but I have no time - needed so it does not mess up the iterations in the first attempt.
    !(because then attmax_special = .true., so we do not set maxatt_iter using the above select case when
    !we go to the next grid point).
    attmax_backup = maxatt_iter

    !The permutation is only done in the box even in the first iteration of VFI and CTS. It might actually
    !make sense to explore everywhere in state space in the first iteration of VFI, not just
    !in the box in the middle - we can afford to there. In particular if we use bspline or some other
    !interpolation scheme which is sensitive to small errors. Otherwise it is probably fine and it saves
    !time. Anyway, in the first stage of CTS we get a rough interpolation only.

    if(pars%VFI_attmax_rnd == 2) then
        !Check if the point is outside of the box. If yes then use probability p2
        !(we should have p2<=p1 but it is not necessary) to determine whether attmax
        !will be set to 1 (we randomise with probability p2 so
        in_the_box = .true.
        if(grids%a1_gr(ind_unfold_all(1,gp_ind)) < pars%VFIperm_a_min(1,1) .or. &
         grids%a1_gr(ind_unfold_all(1,gp_ind)) > pars%VFIperm_a_max(1,1) ) then
            in_the_box = .false.
        elseif(grids%a2_gr(ind_unfold_all(2,gp_ind)) < pars%VFIperm_a_min(1,2) .or. &
         grids%a2_gr(ind_unfold_all(2,gp_ind)) > pars%VFIperm_a_max(1,2) ) then
            in_the_box = .false.
        elseif(grids%rho_gr(ind_unfold_all(3,gp_ind)) < pars%VFIperm_rho_min .or. &
         grids%rho_gr(ind_unfold_all(3,gp_ind)) > pars%VFIperm_rho_max ) then
            in_the_box = .false.
        elseif(grids%aex_gr(ind_unfold_all(4,gp_ind)) < pars%VFIperm_aex_min .or. &
         grids%aex_gr(ind_unfold_all(4,gp_ind)) > pars%VFIperm_aex_max ) then
            in_the_box = .false.
        end if
    end if
    !draw a random number from (0,1) uniformly
    call random_number(VFI_rndnum)
    if(in_the_box) then
        if(VFI_rndnum > pars%VFI_attmax_p1) maxatt_iter = 1
    else
        if(VFI_rndnum > pars%VFI_attmax_p2) maxatt_iter = 1
    end if




    !Check whether PSO should be used at this gridpoint this iteration.
    if(pars%VFI_PSO .and. maxatt_iter > 1) then
        !If we are in the first iteration of VFI (regardless of stage of CTS)
        !If is still used only when maxatt_iter > 1

        !then use it (everywhere where the probability is positive
        call random_number(PSO_rndnum)
        if(PSO_rndnum <= pars%VFI_PSO_p .or. V_ind == 1) then
            !use PSO this iteration
            PSO_this_iter = .true.
            !If we are only using PSO, set maxatt_iter = 1 (this will be restored from backup
            !for the next gridpoint).
            if(pars%VFI_PSO_only) maxatt_iter = 1
        else
            PSO_this_iter = .false.
        end if
    end if

        !Solve the maximization step in VFI at this grid point


        !Just to avoid errors, treat the case of non-allocated x_AAD_img separately.
        if(allocated(x_AAD_img)) then
        call VFI_solve_gp(ind_unfold_all(:,gp_ind),grids,pars,N_a,N_rho,N_aex,pol%x_pol_img(:,gp_ind),&
        V_stub(gp_ind),skip_max,loss_gp,x_AARD_gp,x_MARD_gp,x_AAD_img(gp_ind),maxatt_iter,gp_ind,runtime_limit_unix,PSO_this_iter)
        else
        x_AAD_placeholder = 1.0_dp !this is overwrtten inside so unless we keep setting it to 1, this
        !will be used even if it should not (and what's worse a value corresponding to a different gp
        !policy change will be used).
        call VFI_solve_gp(ind_unfold_all(:,gp_ind),grids,pars,N_a,N_rho,N_aex,pol%x_pol_img(:,gp_ind),&
        V_stub(gp_ind),skip_max,loss_gp,x_AARD_gp,x_MARD_gp,x_AAD_placeholder,maxatt_iter,gp_ind,runtime_limit_unix,PSO_this_iter)
        end if

    !RESTORE attmax from backup.
    maxatt_iter = attmax_backup

        !update total loss at image, and other statistics (for deviation of policy function)
        loss_img_tot = loss_img_tot + loss_gp
        x_MARD_img = max(x_MARD_img,x_MARD_gp)
        x_AARD_img = x_AARD_img + x_AARD_gp !This will be divided by the normalized later. For not it is
        !the addition of the average ARD (across nc_varpar) on the image.
    end do !end cycling over the grid points.

    !After computing its part of value function, every image sends its peace of
    !it to the special image
    call send_V_to_spec(shared_data_V,V_stub,spec_img,min_ind_img,max_ind_img)
    sync all !Sync to make sure all pieces arrived


    sync all !Sync before recovering these numbers at the special image
    if(i_img == spec_img) then
        x_MARD = 0.0_dp
        x_AARD = 0.0_dp
        do img_ind = 1,n_img
            !Maximum absolute (value of) relative difference
            x_MARD = max(x_MARD,x_MARD_img[img_ind])
            !Average absolute (value of) relative difference
            x_AARD = x_AARD + x_AARD_img[img_ind]/(V_num_el)
        end do
    end if

    !Compute average loss across all grid points.
    !loss_img_tot, loss_all_avg
    loss_all_avg = 0.0_dp
    do img_ind = 1,n_img
        loss_all_avg = loss_all_avg + loss_img_tot[img_ind]/(V_num_el)
    end do

    !DEBUG:
    !Compute minimum and maximum value of V_old at this iteration (for debugging)
    if(i_img == spec_img) then

    !if bspline_check_V_convergence = .true. and we are using bspline, then the debug values should be coputed using
    !V_old_copy instead
    if(pars%bspline_check_V_convergence .and. pars%VFI_interpolation_mode == 3) then
        if(.not. (allocated(V_old_copy))) then
        write(*,*) 'Warning: V_old_copy is not allocated when it should be in union.f90.'
        write(*,*) 'V_old_min and V_old_max are set to zero! investigate.'
        debug_val = 0.0_dp
        debug_val2 = 0.0_dp

        else
        debug_val = minval(reshape(V_old_copy,[N_a*N_a*N_rho*N_aex*pars%M_grid]))
        debug_val2 = maxval(reshape(V_old_copy,[N_a*N_a*N_rho*N_aex*pars%M_grid]))


        end if

    else
        debug_val = minval(reshape(V_old,[N_a*N_a*N_rho*N_aex*pars%M_grid]))
        debug_val2 = maxval(reshape(V_old,[N_a*N_a*N_rho*N_aex*pars%M_grid]))
    end if


        write(output_string,*) 'V_old_min = ',debug_val,'V_old_max = ',debug_val2
        call double_output(folder_name,trim(output_string))

        if((debug_val<-1000000000.0_dp .or. debug_val2>1000000000.0_dp) .and. V_ind>1) then
            write(*,*) 'Something bad happened (search for debug_val in code)'
            !error stop
        end if

        write(*,*) 'Test removing V_old_min and V_old_max computation and see how much time we save'

    end if

    !Now image spec_img has all the data. We can't read them directly because of segmentation
    !errors (at least in the version of Intel compiler which I used - 15.0.something,
    !it could also be system-dependent), so
    !on image spec_img we first reshape the data and save the result in a temporary variable to which we allocate
    !memory. Then all other images retrieve data.
    if(i_img == spec_img) then
        if(.not. allocated(shared_data_V%V_tmp)) then
                 allocate(shared_data_V%V_tmp(N_a,N_a,N_rho,N_aex,pars%M_grid))
        end if
        call V_unfold(shared_data_V)
    end if

    !Now that we have the new value function on image spec_img, compute some statistics
    !Maximum absolute relative difference b/w new and old value function
    if(i_img == spec_img) then

    !If we are using bspline, the change in V criterion is computed differently (because V_old was overwritten).
    if(pars%bspline_check_V_convergence .and. pars%VFI_interpolation_mode == 3) then
    V_MARD = maxval(reshape(abs((shared_data_V%V_tmp - V_old_copy)/(abs(V_old)+1.0_dp)),[N_a*N_a*N_rho*N_aex*pars%M_grid]))
    V_AARD = sum(reshape((abs(shared_data_V%V_tmp - V_old_copy)/(abs(V_old)+1.0_dp))/&
        (N_a*N_a*N_rho*N_aex*pars%M_grid),[N_a*N_a*N_rho*N_aex*pars%M_grid]))

    else
    V_MARD = maxval(reshape(abs((shared_data_V%V_tmp - V_old)/(abs(V_old)+1.0_dp)),[N_a*N_a*N_rho*N_aex*pars%M_grid]))
    V_AARD = sum(reshape((abs(shared_data_V%V_tmp - V_old)/(abs(V_old)+1.0_dp))/&
        (N_a*N_a*N_rho*N_aex*pars%M_grid),[N_a*N_a*N_rho*N_aex*pars%M_grid]))
    endif


    !IF debug mode is on, also compute and display additional statistics about value function.
    if(pars%debug_mode) then
        V_avg =  sum(reshape((shared_data_V%V_tmp)/(N_a*N_a*N_rho*N_aex*pars%M_grid),[N_a*N_a*N_rho*N_aex*pars%M_grid]))
    end if

    end if

    !Report some stuff
    if(i_img == spec_img) then

        write(output_string,*) 'V_MaxARD = ',V_MARD
        call double_output(folder_name,trim(output_string))
        write(output_string,*) 'V_AvgARD = ',V_AARD
        call double_output(folder_name,trim(output_string))

        !change in the policy function.
        write(output_string,*) 'x_MaxARD = ',x_MARD,'x_AvgARD = ',x_AARD
        call double_output(folder_name,trim(output_string))

        write(output_string,*) 'loss_all_avg = ',loss_all_avg
        call double_output(folder_name,trim(output_string))

        if(pars%debug_mode) then
        write(output_string,*) 'V_avg = ',V_avg
        call double_output(folder_name,trim(output_string))
        end if
    end if
    
    !Before checking the stopping rule - on image 1 only (so we do not get IO erros when multiple images
    !try to read the same file at the same time) check if the 'stopping file' exists. If yes then
    !set the runtime limit to zero on all images. In VFI it would be sufficient to do this for one image only,
    !but elsewhere (in simulation) it is necessary to change this on all images, so we do it here as well
    !for consistency (no performance implications)
    if(this_image() == 1) then
        call stop_now_check(folder_name,stop_now)
        if(stop_now) then
            runtime_limit_unix = 0.0_dp
            runtime_limit_unix_cts = 0.0_dp
            runtime_limit_unix_how = 0.0_dp
            call double_output(folder_name,'****************************')
            call double_output(folder_name,'File stop_now.txt found. Stopping computations...')
            call double_output(folder_name,'****************************')
        end if
    end if
    syncall
    !all images synchronize the runtime limits with image 1
    if(this_image() /= 1) then
        runtime_limit_unix = runtime_limit_unix[1]
        runtime_limit_unix_cts = runtime_limit_unix_cts[1]
        runtime_limit_unix_how = runtime_limit_unix_how[1]
    end if
    syncall

    !Stopping rule
    !(check only when we're not skipping maximization in Howard algorithm)
    if(i_img==spec_img .and. skip_max == .false.) then
        if(pars%stopping_rule_quantile > 0.0_dp .and. pars%stopping_rule_quantile < 1.0_dp) then
            !Quantile stopping rule
            !Compute the quantile (this can take a while so only do it if quantile stopping rule used)
            allocate(ARD_vector(N_a*N_a*N_rho*N_aex*pars%M_grid))
            ARD_vector=reshape(abs((shared_data_V%V_tmp - V_old)/(abs(V_old)+1.0_dp)),[N_a*N_a*N_rho*N_aex*pars%M_grid])
            ifail = 0
            quant_nag_input(1) = pars%stopping_rule_quantile
            call g01amf(N_a*N_a*N_rho*N_aex*pars%M_grid,ARD_vector,1,quant_nag_input,V_quant_ARD,ifail)
            write(tmp_string1,'(F5.2)') pars%stopping_rule_quantile
            write(tmp_string2,*) V_quant_ARD(1)
            output_string = trim(tmp_string1)//' V_ARD quantile = '//trim(tmp_string2)
            call double_output(folder_name,trim(output_string))

            if(V_quant_ARD(1)<pars%stopping_rule_threshold) then
                VFI_stop_now = .true.
                output_string = 'Stopping rule ('//trim(tmp_string1)//' quant abs rel diff) satisfied, exiting VFI loop.'
                call double_output(folder_name,trim(output_string))
            end if

            if(DA_this_iter) deallocate(ARD_vector)
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

    !If we reached 80% of the runtime (or something else - it is hardcored in file
    !taxhetgr_initialize.f90) do not use Howard's acceleration algorithm anymore.
    if(time_real_acc() > runtime_limit_unix_how .and. pars%crit_how > 0 .and. &
    CTS_ind == pars%CTS_split_times + 1) then

        skip_max = .false. !skip no more maximizations
        !and do not use Howard's algorithm in future iterations
        pars%crit_How = 0
        use_How = .false.
        skip_count_how = 0

        if(i_img == spec_img) then
            call double_output(folder_name,'80% of runtime limit reached. Howard algorithm will not be used anymore.')
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
                crit_value_How = x_MARD
            case(4)
                crit_value_How = x_AARD
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
!on coarse grid). Also don't save if save_iter_multiple is 0, and if sim_only == .true. - this
!would result in using up a lot of HDD space and having to keep deleting value
!functions.
    if((mod(V_ind,pars%save_iter_multiple) == 0 .or. V_ind == pars%VFI_max_iter .or. VFI_stop_now[spec_img])&
        .and. CTS_ind == pars%CTS_split_times + 1 .and. pars%save_iter_multiple/=0 .and. &
        pars%sim_only == .false.) then

        if(i_img == spec_img) then
            call save_V(shared_data_V%V_stub_stacked,grids,N_a,N_rho,N_aex,pars%M_grid,folder_name)
        end if

        !If pars%save_pol, save the optimal choices (policy function).
        if(pars%save_pol) then
            if(i_img == spec_img) then
                write(*,*) 'Saving optimal choices...'
                call runtime_report(initial_time,.true.,elapsed_time_string)
                call double_output(folder_name,trim(elapsed_time_string))
            end if

        !on image spec_img, allocate storage space for guesses from all images
        sync all
        if((i_img == spec_img .or. pars%NDA)) then
            if(.not. allocated(pol%x_pol_all)) allocate(pol%x_pol_all(ncvar_par,V_num_el))
        else !this prevents crashes even though it should not be necessary (specific to compiler)
            if(.not. allocated(pol%x_pol_all)) allocate(pol%x_pol_all(1,1))
        end if
        sync all !after allocating memory

        !Send the individual policy function to image spec_img
        pol[spec_img]%x_pol_all(1:ncvar_par,min_ind_img:max_ind_img) = pol%x_pol_img(1:ncvar_par,min_ind_img:max_ind_img)

        sync all !(after sending)
        if(i_img == spec_img) then
            !Save into file
            call save_policy(pol%x_pol_all,folder_name)
        end if

        sync all !sync before deallocating. Not doing this could lead to some unexpected errors.
        if(DA_this_iter) deallocate(pol%x_pol_all)
            if(i_img == spec_img) then
                write(*,*) 'Optimal choices saved...'
                call runtime_report(initial_time,.true.,elapsed_time_string)
                call double_output(folder_name,trim(elapsed_time_string))
            end if
        end if !(end of branch for saving guess of policy function)
        sync all
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

    call get_V_from_spec(shared_data_V,V_old,spec_img)

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
!_____________end of VFI_______________________________________

!Before the next iteration (or at the end if this is the end indeed), deallocate V_old (on image one, save
!it to a temporary variable if this is not the last iteration)

!If we're not at the last stage of CTS, make sure that all images have a copy of the whole policy function
!from this iteration when we begin the next iteration. This will be used to get an initial guess.
!only do this if pars%CTS_interp == .true. (otherwise we only interpolate over value function and we don't
!need to pass the policy functions).
if(CTS_ind < pars%CTS_split_times + 1 .and. pars%CTS_C_interp) then
    !allocate memory for all of the policy function on all images in the unfolded form, and in the raw form
    !(in which there is one index only for cycling) on image 1.
    if(allocated(pol%x_pol_all)) then
        write(*,*) 'Warning: pol%x_pol_all should not be allocated on any image after end of VFI algorithm.'
        deallocate(pol%x_pol_all)
    end if
    if(allocated(pol%x_pol_all_unf)) then
        write(*,*) 'Warning: pol%x_pol_all_unf should not be allocated on any image after end of VFI algorithm.'
        deallocate(pol%x_pol_all_unf)
    end if
    if(this_image() == spec_img) then
        allocate(pol%x_pol_all(ncvar_par,V_num_el))
    else
        allocate(pol%x_pol_all(1,1))
    end if
    sync all

    !on all images, allocate the unfolded policy function.
    allocate(pol%x_pol_all_unf(ncvar_par,N_a,N_a,N_rho,N_aex,pars%M_grid))
    !Send the policy function part from all images to image 1, so there we get a policy function
    !for all grid points.
    sync all !(make sure the storage space is allocated on all images before sending the policy function)
    !Send the policy function to image 1
    pol[spec_img]%x_pol_all(1:ncvar_par,min_ind_img:max_ind_img) = pol%x_pol_img(1:ncvar_par,min_ind_img:max_ind_img)

    sync all
    !Reshape the policy on image 1 and save it in pol%x_pol_unf (if we try to reshape directly when
    !reading on other images, the program will crash. It's a bug in Intel Fortran compiler (at least it was when
    !I was writing this program, it might have been fixed at some point).
    if(this_image() == spec_img) then
        pol%x_pol_all_unf = reshape(pol%x_pol_all,[ncvar_par,N_a,N_a,N_rho,N_aex,pars%M_grid])
    end if

    sync all
    !Send the policy function to all images (in the unfolded form, this is what we need for interpolation
    if(this_image() /= spec_img) then
        pol%x_pol_all_unf = pol[spec_img]%x_pol_all_unf
    end if
    sync all
    !Now deallocate memory for the policy function in 'folded' form
    deallocate(pol%x_pol_all)
    sync all
    !Note: There are perhaps more sync alls than needed but this should lead to very little time lost (maybe a second)
    !and it could avoid some crashes.

    !Now all images have a copy of the policy function on all grid points in pol%x_pol_all_unf. This is not an issue in terms
    !of memory because it's never from the final stage of CTS algorithm so these are small objects.
end if

!deallocate these before going to the next state of CTS algorithm - the dimensions will be different in the next one.
deallocate(pol%x_pol_img)

if(allocated(x_AAD_img)) deallocate(x_AAD_img) !this is allocated only if rhoendadapt = .true.

!deallocate(optim_wb_img)
deallocate(V_stub)
!If this is the last iteration of CTS algorithm, keep V_old
if(CTS_ind < pars%CTS_split_times + 1) then
    deallocate(V_old)
end if

deallocate(ind_unfold_all)
!This thing could in theory be deallocated in every iteration over V_ind after
!passing data to save memory - but in practice this leads to errors in passing
!data between images on different nodes when the passed arrays are very large
!(order of hundreds of MBs and more)!!! Some bug in CAF implementation...
if(i_img == spec_img .and. allocated(shared_data_V%V_tmp)) deallocate(shared_data_V%V_tmp)

!If this is not the last stage of CTS algorithm, save the value function in vector
!form and the old grids. Do this on image spec_img only.
if(CTS_ind < pars%CTS_split_times + 1) then
    !save the old number of gridpoints
    N_a_tmp = N_a
    N_rho_tmp = N_rho
    N_aex_tmp = N_aex
    if(i_img == spec_img) then
        allocate(V_vect_tmp(N_a_tmp*N_a_tmp*N_rho_tmp*N_aex_tmp*pars%M_grid))
        V_vect_tmp = shared_data_V%V_stub_stacked
    end if
    !grids are needed on all images, not just special image like the value function.

    allocate(grids_tmp%A1_gr(N_a_tmp),grids_tmp%A2_gr(N_a_tmp),grids_tmp%rho_gr(N_rho_tmp),&
    grids_tmp%aex_gr(N_aex_tmp))
    grids_tmp = grids
end if

!Once saved, deallocate the old grids. This sync all should actually not be necessary
!because here, image spec_img only reads its own copy of the data.
sync all
if(allocated(shared_data_V%V_stub_stacked)) deallocate(shared_data_V%V_stub_stacked)

!Keep the grids if this is the last stage of CTS
if(CTS_ind < pars%CTS_split_times + 1) then
    deallocate(grids%A1_gr,grids%A2_gr,grids%rho_gr,grids%aex_gr)
end if
end do !End of CTS algorithm (and thus VFI).

!sync after we got solution of the VFI problem
sync all


!If we are using spline, compute the coefficients here. Alsom write how long it takes.
if(pars%VFI_interpolation_mode == 3) then

sync all !this one should not be necessary?

if (i_img == spec_img) then
    call double_output(folder_name,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    call double_output(folder_name,'Before computing spline coefficients:')
    call runtime_report(initial_time,.true.,elapsed_time_string)
    call double_output(folder_name,trim(elapsed_time_string))
    call double_output(folder_name,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
end if

!initialize the following to 1
inbvx = 1
inbvy = 1
inbvz = 1
inbvq = 1
iloy = 1
iloz = 1
iloq = 1

do m_ind = 1,pars%M_grid
    call db4ink(grids%a1_gr,N_a,grids%a2_gr,N_a,grids%rho_gr,N_rho,grids%aex_gr,N_aex,&
    V_old(:,:,:,:,m_ind),pars%bspline_k_a1,pars%bspline_k_a2,pars%bspline_k_rho,pars%bspline_k_aex,0,&
    t_a1(1:(N_a + pars%bspline_k_a1),m_ind),t_a2(1:(N_a + pars%bspline_k_a2),m_ind),&
    t_rho(1:(N_rho + pars%bspline_k_rho),m_ind),t_aex(1:(N_aex + pars%bspline_k_aex),m_ind),&
    V_old(:,:,:,:,m_ind),bspline_iflag)
end do

if (i_img == spec_img) then
    call double_output(folder_name,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    call double_output(folder_name,'After computing spline coefficients:')
    call runtime_report(initial_time,.true.,elapsed_time_string)
    call double_output(folder_name,trim(elapsed_time_string))
    call double_output(folder_name,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
end if

sync all !just so we don't display output out of order. It doesn't matter because
!all images have to wait for image 1 anyway in simulation, so there will be another sync all
!in a bit.

end if !end of preparing bspline coeffs


!______________________SIMULATION______________________________________________
!The simulation is parallelized in such way that only the first image (CPU) handles all
!the data. The role of other images is only to run maximizations with different
!random initial conditions - this way we can (in a given time) explore more solutions.
!In case of aggregate shocks this is usually not that useful because the algorithms
!tend to converge to a point which is either the same or similar as long as the initial
!conditions are not too far from the solution. But in the case of idiosyncratic shocks,
!this is very useful, as the algoritm often struggles in finding the global optimum
!(given an approximation of value function).

if(this_image() == 1) then
    allocate(CPO_optim_count(num_images()))
    CPO_optim_count = 0
end if


!Test: add a huge constant to the middle point of the grid and see what happens.
!In particular - is it possible to stay in the middle of the grid forever? Deterministic
!case suggests that there are problems with this due to small error in the first period.
if(pars%sim_debug  .and. abs(pars%sim_debug_const) > 1.0E-15_dp) then
    if(this_image() == 1) then
        write(*,*) 'Warning. sim_debug_true and sim_debug_const /= 0.0.'
        write(*,*) 'The results of the simulation are for debugging only, they are not valid!!!'
        writE(*,*) 'index of the midpoint = ',pars%N_a/2 + 1,pars%N_a/2 + 1,pars%N_rho/2 + 1,pars%N_aex/2 + 1
    end if
    !Add the constant to the middle point of the grid. This already assumes that the number of grid points
    !is odd (which should always be the case - otherwise the constant will be added to a point slighly
    !off the middle).
    !(integer division - it should work for the middle of the grid).

    V_old(pars%N_a/2 + 1,pars%N_a/2 + 1,pars%N_rho/2 + 1,pars%N_aex/2 + 1,:) = &
    V_old(pars%N_a/2 + 1,pars%N_a/2 + 1,pars%N_rho/2 + 1,pars%N_aex/2 + 1,:) + pars%sim_debug_const
end if

!

!If we are loading shocks, do so on image 1 only. The other images do not need to know
!the shock realizations as their only role is in parallelization of maximization at a given
!point in state space
if(pars%load_shocks) then
    sync all
    if(this_image() == 1) then
        !allocate memory to the shocks.
        allocate(s_ind_loaded(pars%shocks_T,pars%shocks_N))
        call load_shocks(folder_name,pars,s_ind_loaded)
    end if
    sync all
end if

    !Initialize RNG at some value - if we do it within subroutine that generate Markov chain,
    !etc., there are segmentation faults. Some issue with implementation of CAF by Intel...
    call RANDOM_SEED(SIZE = seed_n)
    allocate(seed_array(seed_n))

    if(pars%rng_seed < 0) then
        call random_seed !if seed negative,just initialize RNG randomly.
    else
        call random_seed(pars%rng_seed) !initialize rng seed for repeatability
    end if
    !call random_seed(pars%rng_seed) !use call random_seed without arguments for random initialization

    !On images other than 1 perform further randomisation based on seed (if we just use call random_seed we might
    !get the same draws on different CPUs - it uses system time, rounded to the closest second I think).
    if(this_image() /= 1) then
        call random_seed
        call random_seed(get = seed_array)
        seed_array = seed_array - i_img*10
        seed_array = max(seed_array,0)
        call random_seed(put = seed_array)
    end if

    if(i_img == spec_img) then
        call double_output(folder_name,'========Computing simulated series========')
    end if





    !What is really important is to allocate memory to all the simulated series only on
    !image 1. Otherwise we may well run out of memory.

    if(this_image() == 1) then
    !Make a copy of initial conditions for assets (so that we don't have to rewrite the
    !following subroutines and can refer to pars%b_init). Then use this copy to
    !generate permutated b_init, and at the end of the simulation algorithm, rewrite
    !pars%b_init by the original content.
    b_init_copy = pars%b_init
    !Also copy the initial external debt of country 1
    b_ex_init_copy = pars%b_ex_init

    !Allocate memory to simulated series
    call sim_result_alloc(sim_result,pars,1)
    !and to the average series
    call sim_result_alloc(sim_avg,pars,1)

    !Allocate memory to the type containing all of the simulations:
    allocate(sim_res_all%SRA(pars%N_sim))

    !Tracking quantiles of state variables (across simulations).
    allocate(sim_res_all%rho_prime_05(pars%T_sim))
    allocate(sim_res_all%rho_prime_95(pars%T_sim))
    allocate(sim_res_all%aex_prime_05(pars%T_sim))
    allocate(sim_res_all%aex_prime_95(pars%T_sim))
    allocate(sim_res_all%a_prime_05(pars%T_sim,I_par))
    allocate(sim_res_all%a_prime_95(pars%T_sim,I_par))



!If perm_eps = 0.0, then we want to simulate the first-period problem only once for every
!possible shock realization. !This saves a lot of time, allows us to get a more accurate solution,
!and reduces fluctuations due to error in first-period problem solution which can sometimes happen.
if(pars%perm_eps == 0.0_dp) then
    !(This is standard, no need to display the message)...
    !write(*,*) 'Perm_eps = 0.0. Computing initial period solution only once per each shock realization'
    do m_ind = 1,M_par
        call solve_FP_problem(V_old,pars%b_init,pars%b_ex_init,m_ind,pars,SS(m_ind),FP_result_allm(m_ind))
    end do

    !Get the guess of choice variables for period t=1
    call gen_guess_t1(FP_result_allm,x_guess_t1)

end if

end if !img == 1


!Get N_sim samples (random Markov chain and permutation of initial conditions).
do sim_ind = 1,pars%N_sim
    if(this_image() == 1) then
    write(*,*) 'Debug in taxhetgr.f90. sim_ind = ',sim_ind
    !every 10th simulation, report the runtime.
    if(mod(sim_ind,10) == 0) then
        if (i_img == spec_img) then
            call runtime_report(initial_time,.true.,elapsed_time_string)
            call double_output(folder_name,trim(elapsed_time_string))
        end if
    end if

    end if

    !Note: I tried the permutations earlier for debugging and understanding the problem. The results
    !reported in the paper are not permutaded (perm_eps = 0).

    !Solve the first-period problem on image 1 only so far (as long as we do not permutate the gains from
    !parallelizing period 0 problem are minimal.

    if(this_image() == 1) then

    !If this is the first simulation, do not permutate. Otherwise do
    if(sim_ind > 1) then
        !Generate permutation (now pars%b_init contains the permutated initial conditions).
        call b_init_permutate(b_init_copy,pars%b_init,b_ex_init_copy,pars%b_ex_init,pars%perm_eps)
    else
        pars%b_init = b_init_copy
        pars%b_ex_init = b_ex_init_copy
    end if

    !s_init_ind_a is the initial index. We either us the loaded one in s_ind_loaded, or we take
    !one randomly or a fixed one  (depending on the parameters).
    if(pars%load_shocks) then
        s_init_ind_a = s_ind_loaded(1,sim_ind)
    else
    !In the 'first period', randomize over initial shock (unless initial shock given). Then pick the appropriate element
    !of FP_result.
    if(pars%s_init_index == -1) then
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

    endif !pars%load_shocks

    !Solve the first-period problem for every possible initial shock (not just the one actually realized).
    !This will be useful for obtaining an initial guess for period t=1 (where we need a state-contingent
    !plan, and using the first-period solution for all shock realizations may be a poor initial guess).
    !(if perm_eps == 0 we already solved this once and since the initial conditions are not permutated,
    !there is no reason to do this again every simulation)
    if(pars%perm_eps > 0.0_dp) then
        if(sim_ind == 1) then
            write(*,*) '****************************'
            write(*,*) 'Using perm_eps > 0.0. In this case it might be worth it parallelizing solve_FP_problem...'
            write(*,*) '****************************'
        end if
        do m_ind = 1,M_par
            call solve_FP_problem(V_old,pars%b_init,pars%b_ex_init,m_ind,pars,SS(m_ind),FP_result_allm(m_ind))
        end do
        !Also get the guess of state-contingent consumption plan for period t=1.
        call gen_guess_t1(FP_result_allm,x_guess_t1)
    end if

    !Now the first-period result is not saved in FP_result but in FP_result_allm(s_init_ind_a)
    !we need to take this into account.

    !Allocate memory to simulated series
    call sim_result_alloc(sim_res_all%SRA(sim_ind),pars,1)

    !pol%x_pol_all_unf = reshape(pol%x_pol_all,[ncvar_par,N_a,N_a,N_rho,N_aex,pars%M_grid])

    end if !(img == 1)

    sync all !sync before simulating series



    !sim_series is called by all images (but image 1 is the only image with all the data).
    !If we are not using interpolation of the policy function then we cannot really reshape x_pol_all
    !(it would result in errors because it is not allocated). Therefore pass a small array otherwise
    !Do the same if we are not on image 1 (when this will usually not be allocated).
    if(pars%sim_interp .and. this_image() == 1) then
        call sim_series(sim_result,FP_result_allm(s_init_ind_a),pars,grids,x_guess_t1,&
        reshape(pol%x_pol_all,[ncvar_par,N_a,N_a,N_rho,N_aex,pars%M_grid]),sim_shared,s_ind_loaded,&
        sim_ind,CPO_optim_count)
    else !images 2,...,n_img always call the subroutine with different arguments
        !These images should never attempt to access the policy function, otherwise a
        !segmentation fault will occurr.
        call sim_series(sim_result,FP_result_allm(s_init_ind_a),pars,grids,x_guess_t1,&
        reshape([1.0_dp],[1,1,1,1,1,1]),sim_shared,s_ind_loaded,sim_ind,CPO_optim_count)
    end if


    if(this_image() == 1) then
    !Save the whole sim_result
    sim_res_all%SRA(sim_ind) = sim_result
    end if

    !If we got here, increase the number of actually completed simulations by 1
    N_sim_done = N_sim_done + 1

    if(this_image() == 1) then
        call stop_now_check(folder_name,stop_now)
        if(stop_now) then
            runtime_limit_unix = 0.0_dp
            runtime_limit_unix_cts = 0.0_dp
            runtime_limit_unix_how = 0.0_dp
            call double_output(folder_name,'****************************')
            call double_output(folder_name,'File stop_now.txt found. Stopping computations...')
            call double_output(folder_name,'****************************')
        end if
    end if
    syncall
    !all images synchronize the runtime limits with image 1
    if(this_image() /= 1) then
        runtime_limit_unix = runtime_limit_unix[1]
        runtime_limit_unix_cts = runtime_limit_unix_cts[1]
        runtime_limit_unix_how = runtime_limit_unix_how[1]
    end if


    !Check runtime. If it is exceeded, do not perform further simulations.
    if(time_real_acc() > runtime_limit_unix) then
        exit
    end if
end do



if(this_image() == 1) then

if(pars%disc_share > 0.0_dp .and. N_sim_done == pars%N_sim) then
!(never discard outliers if we did not complete all simulations that were supposed to be completed).
!In any case this functionality should never be used and I should perhaps remove it in a future update.
if(this_image() == 1) then
    write(*,*) '=================================='
    write(*,*) 'Discarding outliers should not be done in the current version of the program!!!'
    write(*,*) 'Set disc_share = 0.0 in the parameter file.'
    write(*,*) '=================================='
end if

!Discard outliers. A share disc_share of simulations with the most extreme a_prime (in the sense of closeness to bounds)
!will be discarded. The motivation for this is that sometimes numerical optimization
!just yields wrong answers, particularly in the case of multidimensional optimization
!when feasible sets are not convex, and a few bad simulations hitting the bounds can corrupt the results
!quite badly.
!
!This should be used mainly in the case of sampling over initial conditions for asset holdings.
!If we are interested in average effect of shocks, do not use this as it will cause a bias!

!Get number of discarded observations
if(pars%disc_share < 0.0_dp .or. pars%disc_share > 1.0_dp) then
    call double_output(folder_name,'Error: disc_share must be between 0.0 and 1.0! Setting it to 0.')
    pars%disc_share = 0.0_dp
end if
disc_num = floor(pars%N_sim*pars%disc_share)
!Make sure we keep at least one observation
if(disc_num == pars%N_sim) disc_num = pars%N_sim -1

write(tmp_string1,'(I5)') disc_num
write(tmp_string2,*) 'Discarded ', trim(tmp_string1), ' samples with most extreme a_prime.'
call double_output(folder_name,trim(tmp_string2))

call double_output(folder_name,'All saved statistics are computed after discarding!')
call double_output(folder_name,'To obtain the statistics w/o discarding, set disc_share = 0.0 in parameter file.')

!Allocate storage (just for the simulated series, the rest needs to be allocated only at
!the point of calling subroutine sim_stat, and this will not be done for sim_res_all_tmp).
allocate(sim_res_all_tmp%SRA(pars%N_sim - disc_num))

!Copy all of sim_res_all%SRA(sim_ind)%a_slack_all_avg into a single array so we
!can find the series with the minimum value (the worst solution).

allocate(disc_crit(pars%N_sim),disc_indices(disc_num))
do sim_ind = 1,pars%N_sim
    disc_crit(sim_ind) = sim_res_all%SRA(sim_ind)%a_slack_all_avg
end do

call smallest_indices(disc_crit,disc_indices)

!Copy the undiscarded indices into temporary storage
ind_free = 1 !initialize free index in the tmp storage
do sim_ind = 1,pars%N_sim
    !if the index of simulation is not to be discarded copy the sample
    if(.not. any(disc_indices == sim_ind)) then
        sim_res_all_tmp%SRA(ind_free) = sim_res_all%SRA(sim_ind)
        ind_free = ind_free + 1
    end if
end do

!CHANGE pars%N_sim to reflect the number of dropped observations
!This is the easiest way to introduce discarding but care must be taken if
!pars%N_sim is used afterwards. In particular now pars%N_sim is not adjusted
!on other images!!! So far it does not matter but if it's an issue in the future
!introduce pars%N_sim_disc and use that afterwards to avoid confusion.
pars%N_sim = pars%N_sim - disc_num

!Rellocate the part of sim_res_all that depends on pars%N_sim and copy the
!data from temporary storage there
deallocate(sim_res_all%SRA)


allocate(sim_res_all%SRA(pars%N_sim))


!Copy data
sim_res_all%SRA(1:pars%N_sim) = sim_res_all_tmp%SRA(1:pars%N_sim)

!Don't keep the temporary data used for discarding
deallocate(disc_crit,disc_indices)
deallocate(sim_res_all_tmp%SRA)

end if !End of discarding observations.
    !Compute series of average realizations and some other statistics (in the case of M=1 the
    !average series will be the one deterministic series, and the statistics regarding stochasticity of solution
    !will not be interesting).

    if(pars%sim_debug_par) then
        call double_output(folder_name,'~~~~~~SimulationParalelisation debug~~~~~~~~~~')
        call double_output(folder_name,'For every image, display the count of times it found the optimum in sim.')
        call double_output(folder_name,'(With the exception of image 1 this should be roughly uniform)')
        !If not - then some images probably have wrong data such as wrong value function or policy function
        !This could be the case on the supercomputer on different nodes, I need to test this.

        do img_ind = 1,num_images()
            write(output_string,*) 'i_img = ', img_ind, ' opt_count = ',CPO_optim_count(img_ind)
            call double_output(folder_name,trim(output_string))
        end do
        call double_output(folder_name,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    end if
    !CPO_optim_count

    call sim_stat(sim_res_all,pars,sim_avg,N_sim_done)

    !Save statistics of interest. Also this saves simulated series for all samples
    !(for variables of interest - so we can generate fan charts, and we can compute some additional statistics
    !using Matlab).
    call save_stat(sim_res_all,folder_name,pars,N_sim_done)

    !Save the average series realization.
    call save_series(sim_avg,folder_name,pars,N_sim_done)


    call sim_result_alloc(sim_result,pars,0)
    call sim_result_alloc(sim_avg,pars,0)

    deallocate(sim_res_all%SRA)

    deallocate(sim_res_all%rho_prime_05)
    deallocate(sim_res_all%rho_prime_95)
    deallocate(sim_res_all%aex_prime_05)
    deallocate(sim_res_all%aex_prime_95)
    deallocate(sim_res_all%a_prime_05)
    deallocate(sim_res_all%a_prime_95)

pars%b_init = b_init_copy !restore the original value (just in case it is needed again after simulation).
pars%b_ex_init = b_ex_init_copy


end if !End of i_img == 1 (and of simulation)

    if(pars%sim_debug_par) then
    if(this_image() == 1) write(*,*) '----checking that rng seed sum is different----'
    sync all
    call random_seed(get = seed_array)
    write(*,*) 'i_img = ',i_img,'seed_sum = ',sum(seed_array)
    sync all
    if(this_image() == 1) write(*,*) '-------------------------'
    endif


!Display runtime and save it in log file
sync all
if (i_img == spec_img) then
    call runtime_report(initial_time,.true.,elapsed_time_string)
    call double_output(folder_name,trim(elapsed_time_string))
end if

if(allocated(shep_inc)) deallocate(shep_inc)

if(allocated(t_a1)) deallocate(t_a1)
if(allocated(t_a2)) deallocate(t_a2)
if(allocated(t_rho)) deallocate(t_rho)
if(allocated(t_aex)) deallocate(t_aex)


!if(i_img == 1) error stop !sometimes the program doesn't terminate for some reason so stop it manually


1988 end program union
