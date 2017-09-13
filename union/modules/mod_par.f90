module mod_par
!This module contains definition of type par which contains all parameters,
!and of methods working with the type.
!
!This module is not general, i.e., has to be changed for every project depending 
!on what the parameters are.

!We also have to make sure that every time we add/remove a parameter in the type
!par, we add/remove the corresponding reading of that parameter to/from the subroutine
!read_parameters. They should be kept in the same order as in the parameter file
!baseline.txt for clarity but it doesn't make any difference regarding functionality.

use mod_types
use mod_parameters_hard
use mod_IO !general routines used for input/output


implicit none


!For meaning of these parameters, see the parameter file baseline.txt which should contain
!all the parameters, along with comments.
type par
        real(dp) :: runtime_limit
        logical :: debug_mode! = .false.
        logical :: NDA
        logical :: NDA_override !If this is true (default is false) then NDA is not set to be equal to false when sim_only
        logical :: gen_transfers_ss
        logical :: sim_only! = .false.
        logical :: sim_interp
        logical :: sim_takebest
        logical :: sim_PSO
        character(256) :: input_folder! = ''
        logical :: load_initial_guess! = .false.
        integer :: T_sim! = 0
        integer :: N_sim! = 0
        real(dp) :: acc_int
        integer :: acc_samples
        integer :: acc_method
        real(dp) :: perm_eps
        real(dp) :: disc_share
        integer :: RNG_seed! = 0
        real(dp), dimension(1,2) :: mass! = 0.0_dp
        real(dp), dimension(1,2) :: l_max! = 0.0_dp
        real(dp) :: beta! = 0.0_dp
        real(dp), dimension(1,2) :: theta_0! = 0.0_dp
        real(dp), dimension(1,I_par) :: g_share! = 0.0_dp
        logical :: IID_shocks
        real(dp), dimension(M_par,M_par) :: P! = 0.0_dp
        logical :: use_P_sim
        real(dp), dimension(M_par,M_par) :: P_sim
        real(dp), dimension(M_par,I_par) :: G_shock!multiplicative shocks
        real(dp), dimension(M_par,I_par) :: G!actual expenditure realizations (this is not read from parfile but computed after reading paramters).
        real(dp), dimension(M_par,I_par) :: G_copy !A copy of G which will be used if gen_transfers_ss = .true. only.
        integer :: s_init_index! = 0
        logical :: B_min_constant! = .false.
        real(dp), dimension(1,2) :: b_init! = 0.0_dp
        real(dp) :: b_ex_init! = 0.0_dp
        real(dp), dimension(1,2) :: alpha! = 0.0_dp
        real(dp) :: gamma_par, A_par,B_par,sigma_par! = 0.0_dp
        real(dp) :: stopping_rule_quantile, stopping_rule_threshold!=0.0_dp
        integer :: VFI_max_iter! = 0
        integer :: CTS_split_times! = 0
        integer :: CTS_gridpoints_floor! = 0
        integer :: CTS_first_maxatt
        integer :: grids_type
        integer :: N_a, N_rho, N_aex
        logical :: rho_gr_symm
        real(dp) :: rho_min, rho_max! = 0.0_dp
        real(dp), dimension(1,2) :: a_min, a_max! = 0.0_dp
        real(dp) :: aex_min, aex_max! = 0.0_dp
        integer :: VFI_interpolation_mode! = 0

        integer :: bspline_k_a1,bspline_k_a2,bspline_k_rho,bspline_k_aex

        integer :: bspline_N_nonlin
        logical :: bspline_check_V_convergence
        integer :: x_interpolation_mode
        integer :: shep_n
        real(dp) :: shep_floor !the min number that is inverted (distance squared + epsilon is inverted).
        logical :: shep_norm
        integer :: save_iter_multiple! = 10
        real(dp), dimension(1,N_how) :: a_How !row vectors
        integer, dimension(1,N_how) :: b_How
        real(dp) :: c_How
        integer :: crit_How
        logical :: CTS_C_interp! = .false.
        logical :: save_pol! = .false.
        logical :: always_perm_first_iter

        !parameters related to permutations in VFI
        integer :: VFI_attmax
        integer :: VFI_attmax_rnd
        real(dp) :: VFI_attmax_p1
        real(dp) :: VFI_attmax_p2
        real(dp) :: VFI_perm
        !VFI_PSO
        logical :: VFI_PSO,VFI_PSO_only
        real(dp) :: VFI_PSO_p,VFI_PSO_perm
        integer :: VFI_PSO_par

        integer :: sim_PSO_par

        logical :: sim_careful

        real(dp) :: VFIperm_rho_min
        real(dp) :: VFIperm_rho_max
        real(dp), dimension(1,2) :: VFIperm_a_min
        real(dp), dimension(1,2) :: VFIperm_a_max
        real(dp) :: VFIperm_aex_min
        real(dp) :: VFIperm_aex_max


        logical :: LFFC_use_new_guess

        real(dp) :: LFFC_kappa

        !new parameters related to accuracy
        integer :: acc_attmax1, acc_attmax2, acc_maxiter
        integer :: LFFC_attmax, LFFC_maxiter, VFI_opt_maxiter, sim_maxiter, sim_attmax1, sim_attmax2

        integer :: sim_attmax2a,sim_attmax2b !these parameters supercede sim_attmax2 (this can be removed later).
        integer :: opt_subroutine
        real(dp) :: rhoend,e1,e2
        real(dp) :: sim_rhoend,sim_e1,sim_e2
        real(dp) :: sim_perm1,sim_perm2
        real(dp) :: sim_PSO_perm
        logical :: sim_debug, sim_debug_par
        real(dp) :: sim_debug_const
        logical :: sim_FP_quadrilin

        logical :: sim_par_keepattmax !default is false. If this is true then every CPU does the number
        !of maximization attempts prescribed by attmax

        logical :: rhoendadapt

        !The next few parameters are related to loading shocks from file and generating shocks
        logical :: load_shocks,shocks_gen_only
        character(256) :: shocks_path
        integer :: shocks_T, shocks_N

        !parameters related to loss function
        real(dp) :: loss1a,loss2a,loss3a, loss_const

        real(dp) :: LFFC_rnd_a, LFFC_rnd_b, LFFC_rnd_c
        logical :: LFFC_neg_loss

        real(dp), dimension(1,I_par) :: Eg !expected government expenditure in all countries, according
        !to steady-state distribution probabilities. Not read from parameter file directly but computed
        !by read_parameters

        real(dp), dimension(3*M_par - 1) :: x_bu !Upper bounds in maximization over any variables for choice vector x
        !(this is described in definition of type x_guess in mod_union).
        !So far there are no aggregate shocks, etc. and government expenditure is constant, so we only need
        !one bound computed once (unlike in previous projects).
        real(dp), dimension(M_par,I_par) :: max_cons_mat !maximum consumption (in matrix form - used in construction
        !of x_bu).

        !__________________________________________________________________________
        !Some variables which are not loaded from file but computed afterwards - nevertheless
        !it is useful to keep them together with the other parameters.
        real(dp),dimension(M_par,M_par,M_par) :: P_onerow! = 0.0_dp
        !Element (s,:,:) contains row s of the transition matrix P stacked on top of
        !each other M_par times. This is useful for computing conditional expectations efficiently
        !(avoiding loops).

        integer :: M_grid !the actual number of grid points for shock realization used in constructing
        !the grid. if IID_shocks = .true., then this will be 1, otherwise it will be M_par.

        !Constraints in matrix form - We compute them once so we don't have to stack them on top of each other
        !every time they are used. This allows us to avoid do cycles, and saves a bit of time (these constraints
        !are used very deep in the program, in subroutine check_constr).
        real(dp), dimension(M_par,I_par) :: l_max_mat,a_min_mat,a_max_mat
        real(dp), dimension(M_par,1) :: rho_min_vec,rho_max_vec,aex_min_vec,aex_max_vec

        real(dp), dimension(4) :: shep_norm_val !values used in Shepard's interpolation normalization if shep_norm == .true.
        !(not read from parameter file but computed by subroutine read_parameters)

end type par

contains

subroutine read_parameters(parameter_file,parameters)
        type(par), intent(out) :: parameters
        character(len=*), intent(in) :: parameter_file
        
        character(256) :: file_path !will contain the whole path to the parameter file.
        
        integer :: i,j !index for a do loop

        real(dp), dimension(1,I_par) :: g_avg !average government expenditure.
        real(dp) :: p1 !steady-state probability of being in state with index 1 (if M_par = 2)

        !file_path: we assume that the parameter files are in subfolder parameters.
        file_path = 'parameters/'//trim(parameter_file)
        
        call GetReal('runtime_limit',file_path,parameters%runtime_limit,11.5_dp)
        call GetLogical('debug_mode',file_path,parameters%debug_mode,.false.)
        call GetLogical('NDA',file_path,parameters%NDA,.false.)
        call GetLogical('NDA_override',file_path,parameters%NDA_override,.false.)

        call GetLogical('gen_transfers_ss',file_path,parameters%gen_transfers_ss,.false.)
        call GetLogical('sim_only',file_path,parameters%sim_only,.false.)
        call GetLogical('sim_interp',file_path,parameters%sim_interp,.false.)
        call GetLogical('sim_takebest',file_path,parameters%sim_takebest,.true.)

        call GetLogical('load_shocks',file_path,parameters%load_shocks,.false.)
        call GetLogical('shocks_gen_only',file_path,parameters%shocks_gen_only,.false.)
        !Only load this if load_shocks == .true. (default values for strings were not implemented)
        if(parameters%load_shocks) then
            call GetString('shocks_path',file_path,parameters%shocks_path)
        else
            parameters%shocks_path = ''
        end if
        call GetInt('shocks_N',file_path,parameters%shocks_N,10)
        call GetInt('shocks_T',file_path,parameters%shocks_T,10)

        call GetLogical('sim_PSO',file_path,parameters%sim_PSO,.false.)

        call GetString('input_folder',file_path,parameters%input_folder)
        call GetLogical('load_initial_guess',file_path,parameters%load_initial_guess,.false.)
        call GetInt('T_sim',file_path,parameters%T_sim,10)
        call GetInt('N_sim',file_path,parameters%N_sim,10)
        call GetReal('acc_int',file_path,parameters%acc_int,0.1_dp)
        call GetInt('acc_samples',file_path,parameters%acc_samples,10)
        call GetInt('acc_method',file_path,parameters%acc_method,1)
        call GetReal('perm_eps',file_path,parameters%perm_eps,0.01_dp)
        call GetReal('disc_share',file_path,parameters%disc_share,0.05_dp)
        call GetInt('RNG_seed',file_path,parameters%RNG_seed,-1)
        call GetRealMat('mass',file_path,parameters%mass)
        call GetRealMat('l_max',file_path,parameters%l_max)
        call GetReal('beta',file_path,parameters%beta,0.95_dp)
        call GetRealMat('theta_0',file_path,parameters%theta_0)
        call GetRealMat('g_share',file_path,parameters%g_share)
        call GetLogical('IID_shocks',file_path,parameters%IID_shocks,.false.)
        call GetRealMat('P',file_path,parameters%P)
        call GetLogical('use_P_sim',file_path,parameters%use_P_sim,.false.)
        call GetLogical('always_perm_first_iter',file_path,parameters%always_perm_first_iter,.false.)

        call GetLogical('sim_par_keepattmax',file_path,parameters%sim_par_keepattmax,.false.)

        call GetInt('sim_PSO_par',file_path,parameters%sim_PSO_par,50)

        !Only load P_sim if use_P_sim = true, otherwise set it equal to P
        if(parameters%use_P_sim) then
            call GetRealMat('P_sim',file_path,parameters%P_sim)
        else
            parameters%P_sim = parameters%P
        end if

        call GetRealMat('G_shock',file_path,parameters%G_shock)
        call GetInt('s_init_index',file_path,parameters%s_init_index,1)
        call GetLogical('LFFC_neg_loss',file_path,parameters%LFFC_neg_loss,.true.)
        !The default option is false (it is not quite debugged at the moment, the
        !old algorithm seems to perform better despite its simplicity).
        call GetLogical('LFFC_use_new_guess',file_path,parameters%LFFC_use_new_guess,.false.)

        call GetLogical('sim_careful',file_path,parameters%sim_careful,.false.)


        call GetRealMat('b_init',file_path,parameters%b_init)
        call GetReal('b_ex_init',file_path,parameters%b_ex_init,0.00_dp)
        call GetRealMat('alpha',file_path,parameters%alpha)
        call GetReal('gamma_par',file_path,parameters%gamma_par,2.0_dp)
        call GetReal('A_par',file_path,parameters%A_par,1.0_dp)
        call GetReal('B_par',file_path,parameters%B_par,1.5_dp)
        call GetReal('sigma_par',file_path,parameters%sigma_par,1.0_dp)
        call GetReal('stopping_rule_quantile',file_path,parameters%stopping_rule_quantile,1.0_dp)
        call GetReal('stopping_rule_threshold',file_path,parameters%stopping_rule_threshold,0.0000001_dp)
        call GetInt('VFI_max_iter',file_path,parameters%VFI_max_iter,1)
        call GetInt('CTS_split_times',file_path,parameters%CTS_split_times,0)
        call GetInt('CTS_gridpoints_floor',file_path,parameters%CTS_gridpoints_floor,5)
        call GetInt('CTS_first_maxatt',file_path,parameters%CTS_first_maxatt,50)
        call GetInt('grids_type',file_path,parameters%grids_type,1)
        call GetInt('N_a',file_path,parameters%N_a,10)
        call GetInt('N_rho',file_path,parameters%N_rho,10)
        call GetInt('N_aex',file_path,parameters%N_aex,10)
        call GetLogical('rho_gr_symm',file_path,parameters%rho_gr_symm,.false.)
        call GetReal('rho_min',file_path,parameters%rho_min,0.1_dp)
        call GetReal('rho_max',file_path,parameters%rho_max,2.0_dp)
        call GetRealMat('a_min',file_path,parameters%a_min)
        call GetRealMat('a_max',file_path,parameters%a_max)
        call GetReal('aex_min',file_path,parameters%aex_min,-5.0_dp)
        call GetReal('aex_max',file_path,parameters%aex_max,5.0_dp)
        call GetInt('VFI_interpolation_mode',file_path,parameters%VFI_interpolation_mode,1)


        call GetInt('bspline_k_a1',file_path,parameters%bspline_k_a1,2)
        call GetInt('bspline_k_a2',file_path,parameters%bspline_k_a2,2)
        call GetInt('bspline_k_rho',file_path,parameters%bspline_k_rho,2)
        call GetInt('bspline_k_aex',file_path,parameters%bspline_k_aex,2)

        call GetInt('bspline_N_nonlin',file_path,parameters%bspline_N_nonlin,5)

        call GetLogical('bspline_check_V_convergence',file_path,parameters%bspline_check_V_convergence,.true.)

        call GetInt('x_interpolation_mode',file_path,parameters%x_interpolation_mode,-1)
        !default value of x_interpolation_mode is -1 (in this case the value given in
        !VFI_interpolation_mode is used instead).



        call GetInt('shep_n',file_path,parameters%shep_n,4)
        call GetReal('shep_floor',file_path,parameters%shep_floor,1.0E-30_dp)

        call GetLogical('shep_norm',file_path,parameters%shep_norm,.false.)


        call GetInt('save_iter_multiple',file_path,parameters%save_iter_multiple,1000)
        call GetInt('crit_How',file_path,parameters%crit_How,10)
        call GetRealMat('a_How',file_path,parameters%a_How)
        call GetIntMat('b_How',file_path,parameters%b_How)
        call GetReal('c_How',file_path,parameters%c_How,0.1_dp)
        call GetLogical('CTS_C_interp',file_path,parameters%CTS_C_interp,.true.)
        call GetLogical('save_pol',file_path,parameters%save_pol,.false.)
        call GetInt('acc_attmax1',file_path,parameters%acc_attmax1,100)
        call GetInt('acc_attmax2',file_path,parameters%acc_attmax2,10)
        call GetInt('acc_maxiter',file_path,parameters%acc_maxiter,500)
        call GetInt('LFFC_attmax',file_path,parameters%LFFC_attmax,10)
        call GetInt('LFFC_maxiter',file_path,parameters%LFFC_maxiter,200)
        call GetInt('VFI_opt_maxiter',file_path,parameters%VFI_opt_maxiter,200)
        call GetInt('sim_maxiter',file_path,parameters%sim_maxiter,200)
        call GetInt('sim_attmax1',file_path,parameters%sim_attmax1,100)
        !sim_attmax2 should be redundant
        call GetInt('sim_attmax2',file_path,parameters%sim_attmax2,10)
        !If these are left at the default values we can expect very poor per
        call GetInt('sim_attmax2a',file_path,parameters%sim_attmax2a,0)
        call GetInt('sim_attmax2b',file_path,parameters%sim_attmax2b,1)
        !The default values of the following parameters are important because these parameters
        !were added recently so they are not in many of the older parameter files.
        call GetInt('VFI_attmax',file_path,parameters%VFI_attmax,1)
        call GetInt('VFI_attmax_rnd',file_path,parameters%VFI_attmax_rnd,0)
        call GetReal('VFI_attmax_p1',file_path,parameters%VFI_attmax_p1,1.0_dp)
        !default prob of randomisation outside of the box will be 0.
        call GetReal('VFI_attmax_p2',file_path,parameters%VFI_attmax_p2,0.0_dp)
        call GetReal('VFI_perm',file_path,parameters%VFI_perm,0.05_dp)

        call GetLogical('VFI_PSO',file_path,parameters%VFI_PSO,.false.)
        call GetLogical('VFI_PSO_only',file_path,parameters%VFI_PSO_only,.true.)
        call GetReal('VFI_PSO_p',file_path,parameters%VFI_PSO_p,1.0_dp)
        call GetInt('VFI_PSO_par',file_path,parameters%VFI_PSO_par,50)
        call GetReal('VFI_PSO_perm',file_path,parameters%VFI_PSO_perm,0.1_dp)


        !Only read the following if VFI_attmax_rnd > 0 (otherwise they are not supposed
        !to be used and most likely are not present in old parameter files - this would lead to crashes)
        if(parameters%VFI_attmax_rnd > 0) then

            call GetReal('VFIperm_rho_min',file_path,parameters%VFIperm_rho_min,0.95_dp)
            call GetReal('VFIperm_rho_max',file_path,parameters%VFIperm_rho_max,1.05_dp)
            call GetReal('VFIperm_aex_min',file_path,parameters%VFIperm_aex_min,-1.0_dp)
            call GetReal('VFIperm_aex_max',file_path,parameters%VFIperm_aex_max,1.0_dp)
            call GetRealMat('VFIperm_a_min',file_path,parameters%VFIperm_a_min)
            call GetRealMat('VFIperm_a_max',file_path,parameters%VFIperm_a_max)

        end if

        call GetReal('LFFC_kappa',file_path,parameters%LFFC_kappa,0.5_dp)

        call GetReal('sim_perm1',file_path,parameters%sim_perm1,0.25_dp)
        call GetReal('sim_perm2',file_path,parameters%sim_perm2,0.25_dp)
        call GetReal('sim_PSO_perm',file_path,parameters%sim_PSO_perm,0.2_dp)
        call GetLogical('sim_FP_quadrilin',file_path,parameters%sim_FP_quadrilin,.false.)

        call GetLogical('sim_debug',file_path,parameters%sim_debug,.false.)

        call GetLogical('sim_debug_par',file_path,parameters%sim_debug_par,.true.)
        call GetReal('sim_debug_const',file_path,parameters%sim_debug_const,0.0_dp)

        call GetInt('opt_subroutine',file_path,parameters%opt_subroutine,2)
        call GetReal('rhoend',file_path,parameters%rhoend,0.000001_dp)

        !Start with a default value = .false., but if it works (which I think it will,
        !change it to .true. at some point).
        call GetLogical('rhoendadapt',file_path,parameters%rhoendadapt,.false.)

        call GetReal('e1',file_path,parameters%e1,0.000000001_dp)
        call GetReal('e2',file_path,parameters%e2,0.0000000001_dp)
        call GetReal('sim_rhoend',file_path,parameters%sim_rhoend,0.000001_dp)
        call GetReal('sim_e1',file_path,parameters%sim_e1,0.000000001_dp)
        call GetReal('sim_e2',file_path,parameters%sim_e2,0.0000000001_dp)

        call GetReal('loss1a',file_path,parameters%loss1a,100.0_dp)
        call GetReal('loss2a',file_path,parameters%loss2a,1000.0_dp)
        call GetReal('loss3a',file_path,parameters%loss3a,10000.0_dp)
        call GetReal('loss_const',file_path,parameters%loss_const,0.0_dp)


        call GetReal('LFFC_rnd_a',file_path,parameters%LFFC_rnd_a,0.1_dp)
        call GetReal('LFFC_rnd_b',file_path,parameters%LFFC_rnd_b,0.1_dp)
        call GetReal('LFFC_rnd_c',file_path,parameters%LFFC_rnd_c,0.1_dp)

        !IF x_interpolation_mode == -1, then use the value of VFI_interpolation_mode instead,
        !unless it is 3, in which case use 1 (bspline not yet implemented for policy function
        !as it's probably not worth it but I might to implement it later for simulations)
        if(parameters%x_interpolation_mode == -1) then
            if(parameters%VFI_interpolation_mode /= 3) then
                parameters%x_interpolation_mode = parameters%VFI_interpolation_mode
            else
                parameters%x_interpolation_mode = 1
            end if

        end if

        !pre-compute a few things which are used very often so are conveniently stored in the parameters object

        !If shocks are IID, rewrite all rows of the transition matrix with the first row. This
        !makes the code more efficient that having to check ever time whether shocks are iid when computing
        !conditional expectations, or having to use min function.
        if(parameters%iid_shocks .and. M_par > 1) then
            do i=2,M_par
                parameters%P(i,:) = parameters%P(1,:)
            end do
        end if

        !the following is useful for computing expectations using linear algebra (fast)
        do i=1,M_par !index of the row that we stack
            do j=1,M_par !stack it M times on top of each other
                parameters%P_onerow(i,j,:) = parameters%P(i,:)
            end do
        end do

        !If we are using symmetry for rho grid, then set rho_max = 1/rho_min. This is just a precaution
        !to avoid errors elsewhere (so that the end point of the grid corresponds to the parameter rho_max
        !in type pars).
        if(parameters%rho_gr_symm) then
            parameters%rho_max = 1.0_dp/parameters%rho_min
        end if

        !Compute the all the government shock realizations.
        !cycle over countries
        g_avg = parameters%g_share * parameters%l_max * parameters%theta_0
        !(G contains the actual government expenditure in both countries for all shock realizations).
        do i=1,M_par
            parameters%G(i,:) = g_avg(1,:) * parameters%G_shock(i,:)
        end do

        !Also compute expected government expenditures in both countries in the long run (steady-state distribution).
        !Call this Eg
        !First get the steady-state probability (assumes M_par = 2)
        if(M_par > 2) then
            if(this_image()==1) then
                write(*,*) 'Error: subroutine read_parameters must be generalized for M_par > 2.'
            end if
            sync all
            error stop
        end if
        if(M_par == 2) then
            p1 = parameters%P(2,1)/(1.0_dp + parameters%P(2,1) - parameters%P(1,1))
            parameters%Eg = p1 * parameters%G(1:1,1:I_par) + (1.0_dp-p1) * parameters%G(2:2,1:I_par)
        else !M_par = 1
            p1 = 1.0_dp
            parameters%Eg = parameters%G(1:1,1:I_par)
        end if


        !IF we are interested in generating a steady state of a model with transfers,
        !then overwrite the government expenditures in each country by the overall expenditures
        !divided between the countries equally in each state. If the size of the coountries is
        !different, we would get a total level of expenditures (by multiplying the individual
        !expenditures by masses), then divide it in such way that the per capita expenditures
        !in each country are the same in every state. This can be done later, it is a low
        !priority functionality. The same goes for other generalizations such as non-identical
        !countries (maybe then we would not want the per capita expenditures the same
        !and we would instead have to solve for a SS in a big country - but this is the easiest
        !way to do it using the currently written code).
        !First back up the current government expenditure matrix
        if(parameters%gen_transfers_ss) then
        parameters%G_copy = parameters%G

        !For each state, get the average expenditure and us it as expenditure in both countries
        !(this assumes I_par = 2)
        parameters%G(1,:) = (parameters%G(1,:) + parameters%G(2,:))/2.0_dp
        parameters%G(2,:) = parameters%G(1,:)


        !Check that in every state the government expenditure (which is now the same in all countries)
        !is the same
        do j = 2,M_par
            if(parameters%G(j,1) /= parameters%G(1,1)) then
                error stop('Total government expenditures must be constant across states if get_transfers_ss = .true.')
            end if
        end do

        !Set this true so that we go directly to simulation.
        parameters%sim_only = .true.

        !Also in this case do not use anything that would result in overwriting the initial guess in
        !subroutine sim_series.
        parameters%sim_interp = .false.
        parameters%sim_pso = .false.

        end if

        !Compute bounds for parameters
        call prep_bounds(parameters%l_max,parameters%a_min,parameters%a_max,parameters%rho_min,parameters%rho_max&
        ,parameters%aex_min,parameters%aex_max,parameters%l_max_mat,parameters%a_min_mat&
        ,parameters%a_max_mat,parameters%rho_min_vec,parameters%rho_max_vec,parameters%aex_min_vec,parameters%aex_max_vec)

        !Compute values for normalization in Shepard's interpolation
        !So far this assumes equispaced grid and also it assumes that in CTS, the
        !change in the number of grid points is proportional (in general, in the first
        !stage of CTS this might not be a great approximation if the CTS_floors
        !are used - I should check that later). But so far this is just for experimenting
        !in simulations and second stage of CTS (or no CTS at all so it should
        !be fine). Also if shep_norm and grids are not equispaced - display an error message
        if(parameters%shep_norm .and. parameters%grids_type /= 1) then
            if(this_image()==1) then
            write(*,*) 'ERROR: we cannot use shep_norm when grids are not equispaced. Generalization needed.'
            write(*,*) 'See file mod_par.f90'
            !Basically the generalization would be that the normalization is the distance
            !b/w grid points which changes depending on where in the grid we are.
            error stop
            end if

        end if
        parameters%shep_norm_val(1) = (parameters%a_max(1,1) - parameters%a_min(1,1))&
        /real(parameters%N_a,dp)!a1
        parameters%shep_norm_val(2) = (parameters%a_max(1,2) - parameters%a_min(1,2))&
        /real(parameters%N_a,dp) !a2
        parameters%shep_norm_val(3) = (parameters%rho_max - parameters%rho_min)&
        /real(parameters%N_rho,dp) !rho
        parameters%shep_norm_val(4) = (parameters%aex_max - parameters%aex_min)&
        /real(parameters%N_aex,dp) !a_ex


        !(reminder - NDA is the parameters which prevents problems when loading policy
        !function on different nodes on HPC).
        !if sim_only, then set NDA = .false. In this case only image 1 needs the value
        !of policy functions so it should not create issues.
        !and we need to keep the memory unoccupied. Introduce parameter NDA_override
        !(so that I can run simulation on multiple nodes on HPC).
        !Even so, there might now be some issues on different nodes due to asymmetric allocation of policy function
        !(even though it is not used). So perhaps it is safer to use sim_only on 1 node only
        !and in case of multiple nodes just use NDA_override.
        if(parameters%sim_only .and. (.not.  parameters%NDA_override)) then
            parameters%NDA = .false.
        end if

end subroutine

!Subroutine prep_bounds converts some of the bounds which are not dependent
!on state into state-by-state matrices of bounds. These are then used for checking
!constraints more efficiently than using loops or if clauses.
subroutine prep_bounds(l_max,a_min,a_max,rho_min,rho_max,aex_min,aex_max,l_max_mat,&
a_min_mat,a_max_mat,rho_min_vec,rho_max_vec,aex_min_vec,aex_max_vec)
    real(dp), dimension(1,I_par), intent(in) :: l_max,a_min,a_max
    real(dp), intent(in) :: rho_min,rho_max,aex_min,aex_max

    real(dp), dimension(M_par,I_par), intent(out) :: l_max_mat,a_min_mat,a_max_mat
    real(dp), dimension(M_par,1), intent(out) :: rho_min_vec,rho_max_vec, aex_min_vec,aex_max_vec

    integer :: i

    do i=1,I_par
        l_max_mat(:,i) = l_max(1,i)
    end do
    do i=1,I_par
        a_min_mat(:,i) = a_min(1,i)
    end do
    do i=1,I_par
        a_max_mat(:,i) = a_max(1,i)
    end do

    rho_min_vec(:,1) = rho_min
    rho_max_vec(:,1) = rho_max

    aex_min_vec(:,1) = aex_min
    aex_max_vec(:,1) = aex_max


end subroutine prep_bounds

!Subroutine check_parameters checks certain conditions on parameters and if these are
!violated it either diplays a warning or an error message and stops execution of the program
!(depending on severity of the problem).

subroutine check_parameters(parameters,folder_name)
        type(par), intent(in) :: parameters
        character(len=*), intent(in) :: folder_name !need this so we know where the log file is
        integer :: i
        real(dp) :: sum_check
        character(256) :: error_message
        logical :: critical_error_found
        
        critical_error_found = .false. !If this is true, the program will stop execution
        
        !Check whether the mass vector is non-negative and sums up to 1
        sum_check = 0.0_dp
        do i=1,I_par
                if(parameters%mass(1,i) < 0.0_dp) then
                        write(error_message,*) 'Mass of agent ',i,' is less than one.'
                        call double_output(folder_name, error_message)
                        critical_error_found = .true.
                end if
                sum_check = sum_check + parameters%mass(1,i)
        end do
        if(sum_check /= 1.0_dp) then
                write(error_message,*) 'Mass does not sum up to one.'
                call double_output(folder_name, error_message)
                critical_error_found = .true.
        end if

        !If the grid for rho is symmetric around 1, we need rho_min in (0,1).
        if(parameters%rho_gr_symm .and. (parameters%rho_min <= 0.0_dp .or. parameters%rho_min >= 1.0_dp)) then
                write(error_message,*) 'When rho_gr_symm = .true. in parameter file, rho_min must be in (0,1) interval.'
                call double_output(folder_name, error_message)
                critical_error_found = .true.
        end if


        !Lower bound for grid for rho can't be negative or zero!
        if(parameters%rho_min <= 0.0_dp) then
                write(error_message,*) 'lower bound for rho grid must be greater than zero.'
                call double_output(folder_name, error_message)
                critical_error_found = .true.
        end if
        
        !The number of points per dimension used in Shepard's interpolation cannot be greater than the
        !number of grid point per dimension globally.
        if(parameters%shep_n > parameters%N_a .or. parameters%shep_n > parameters%N_rho .or. &
            parameters%shep_n > parameters%N_aex .or. &
            parameters%shep_n > parameters%CTS_gridpoints_floor .and. parameters%CTS_split_times > 0) then
                write(error_message,*) 'Error: shep_n cannot be greater than N_a, N_rho, N_aex, or CTS_gridpoints_floor'
                call double_output(folder_name, error_message)
                critical_error_found = .true.
        end if

        !Also, shep_n as to be greater than or equal to 2
        if(parameters%shep_n < 2) then
                write(error_message,*) 'Error: parameter shep_n has to be greater than or equal to 2.'
                call double_output(folder_name, error_message)
                critical_error_found = .true.
        end if


        if(parameters%sim_perm1 < 0.0_dp .or. parameters%sim_perm1 >= 1.0_dp) then
                write(error_message,*) 'Error: parameter sim_perm1 has to be in [0,1) interval.'
                call double_output(folder_name, error_message)
                critical_error_found = .true.
        end if
        if(parameters%sim_perm2 < 0.0_dp .or. parameters%sim_perm2 >= 1.0_dp) then
                write(error_message,*) 'Error: parameter sim_perm2 has to be in [0,1) interval.'
                call double_output(folder_name, error_message)
                critical_error_found = .true.
        end if


        if(parameters%VFI_attmax_rnd < 0 .or. parameters%VFI_attmax_rnd > 2) then
                write(error_message,*) 'Error: Parameter VFI_attmax_rnd must be one of 0,1,2.'
                call double_output(folder_name, error_message)
                critical_error_found = .true.
        end if


        if(abs(parameters%LFFC_rnd_b) > 0.3 .and. M_par /= 2) then
                write(error_message,*) 'For M_par not equal to 2, LFFC_rnd_b should be at most 0.3 in absolute value.'
                call double_output(folder_name, error_message)
                write(error_message,*) '(LFFC should be adjusted at some point to generalize the M_par 2 case)'
                call double_output(folder_name, error_message)
        end if


        if(parameters%load_initial_guess .and. num_images() > 16 .and. (.not. parameters%NDA)) then
                call double_output(folder_name,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
                write(error_message,*) 'Warning: NDA should be true in parameter file when load_initial_guess and the number of images is greater than 16'
                !(the 16 images is a proxy for running the program on more than 1 node)
                !This is not a critical error - it might just lead to images on other node than 1 to have the
                !wrong value of policy function (but value function should be fine).
                call double_output(folder_name, error_message)
                call double_output(folder_name,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        end if

        !Related to loading shocks. If N_sim or T_sim is greater than shocks_N or shocks_T, then
        !display a warning and ignore shocks.
        if(parameters%load_shocks) then
        if(parameters%T_sim > parameters%shocks_T .or. parameters%N_sim > parameters%shocks_N) then
                write(error_message,*) 'When loading shocks, we need shocks_T >= T_sim, and shocks_N >= N_sim'
                call double_output(folder_name, error_message)
                write(error_message,*) 'Setting load_shocks = .false. (ignoring shocks file)'
                call double_output(folder_name, error_message)
        end if

        end if


!                logical :: load_shocks,shocks_gen_only
!        character(256) :: shocks_path
!        integer :: shocks_T, shocks_N


        if(critical_error_found) then
                call double_output(folder_name,'Critical error found when checking parameters...')
                error stop
        end if
        
end subroutine check_parameters

end module mod_par
