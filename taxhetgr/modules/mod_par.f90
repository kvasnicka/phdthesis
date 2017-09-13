module mod_par
!This module contains definition of type par which contains all parameters,
!and of methods working with the type.


use mod_types
use mod_parameters_hard
use mod_IO !general routines used for input/output


implicit none


!Type par assumes 2 agents but possibly stochastic environment.
!
!In this file, there are not many comments - more extensive comments are
!in sample parameter file baseline.txt.
type par
        real(dp) :: runtime_limit
        logical :: search_for_bounds! = .false.
        logical :: debug_mode! = .false.
        logical :: static_environment! = .false.
        logical :: sim_only! = .false.
        character(256) :: input_folder! = ''
        logical :: load_initial_guess! = .false.
        integer :: T_sim! = 0
        integer :: N_sim! = 0

        real(dp) :: acc_int
        integer :: acc_samples
        integer :: acc_method

        real(dp) :: perm_eps
        real(dp) :: disc_share
        integer :: T_sim_max! = 0
        integer :: RNG_seed! = 0
        real(dp), dimension(1,2) :: mass! = 0.0_dp
        real(dp), dimension(1,2) :: l_max! = 0.0_dp
        real(dp) :: beta! = 0.0_dp
        real(dp), dimension(1,2) :: theta_0! = 0.0_dp
        real(dp), dimension(1,2) :: xi_theta! = 0.0_dp
        real(dp) :: k_gov! = 0.0_dp
        real(dp), dimension(1,M_par) :: S! = 0.0_dp
        real(dp), dimension(M_par,M_par) :: P! = 0.0_dp
        integer :: s_init_index! = 0
        real(dp), dimension(1,2) :: b_min, b_max! = 0.0_dp
        real(dp) :: B_min_kappa! = 0.0_dp
        real(dp) :: B_min_number! = 0.0_dp
        logical :: B_min_constant! = .false.
        real(dp), dimension(1,2) :: b_init! = 0.0_dp
        real(dp), dimension(1,2) :: alpha! = 0.0_dp
        real(dp) :: gamma_par, A_par,B_par! = 0.0_dp
        logical :: interior_solution_only! = .false.
        logical :: one_fc_suff! = .false.
        real(dp) :: stopping_rule_quantile, stopping_rule_threshold!=0.0_dp
        integer :: VFI_max_iter! = 0
        integer :: VFI_glob_opt_iter! = 0
        logical :: assume_interior_optimality! = .false.
        integer :: LFFC_max_ind! = 10
        real(dp) :: LFFC_min_share! = 0.001_dp
        real(dp) :: LFFC_max_share! = 0.999_dp
        integer :: CTS_split_times! = 0
        integer :: CTS_gridpoints_floor! = 0
        integer :: N_a, N_rho, N_T
        real(dp) :: rho_min, rho_max! = 0.0_dp
        real(dp), dimension(1,2) :: a_min, a_max! = 0.0_dp
        integer :: VFI_interpolation_mode! = 0
        real(dp) :: rho_search_mid, rho_search_max_diff, rho_search_floor! = 0.0_dp
        real(dp), dimension(1,2) :: a_search_min, a_search_max! = 0.0_dp
        integer :: N_search_a, N_search_rho! = 0
        logical :: check_corners_only! = .false.
        logical :: stop_when_found! = .false.
        integer :: save_iter_multiple! = 10
        real(dp), dimension(1,N_how) :: a_How !row vectors
        integer, dimension(1,N_how) :: b_How
        real(dp) :: c_How
        integer :: crit_How
        logical :: CTS_C_interp! = .false.
        logical :: save_C_guess! = .false.

        !new parameters related to accuracy
        integer :: acc_attmax1, acc_attmax2, acc_maxiter
        integer :: LFFC_attmax, LFFC_maxiter, VFI_opt_maxiter, sim_maxiter, sim_attmax1, sim_attmax2
        integer :: opt_subroutine
        real(dp) :: rhoend,e1,e2
        real(dp) :: sim_rhoend,sim_e1,sim_e2

        !__________________________________________________________________________
        !Some variables which are not loaded from file but computed afterwards - nevertheless
        !it is useful to keep them together with the other parameters.
        real(dp),dimension(M_par,M_par,M_par) :: P_onerow! = 0.0_dp
        !Element (s,:,:) contains row s of the transition matrix P stacked on top of
        !each other M_par times. This is useful for computing conditional expectations efficiently
        !(avoiding loops).

        real(dp) :: rel_growth! = 0.95_dp
end type par

contains

subroutine read_parameters(parameter_file,parameters)
        type(par), intent(out) :: parameters
        character(len=*), intent(in) :: parameter_file
        
        character(256) :: file_path !will contain the whole path to the parameter file.
        
        integer :: i,j !index for a do loop

        !file_path: we assume that the parameter files are in subfolder parameters.
        file_path = 'parameters/'//trim(parameter_file)
        
        call GetReal('runtime_limit',file_path,parameters%runtime_limit,11.5_dp)
        call GetLogical('search_for_bounds',file_path,parameters%search_for_bounds)
        call GetLogical('debug_mode',file_path,parameters%debug_mode)
        call GetLogical('static_environment',file_path,parameters%static_environment)
        call GetLogical('sim_only',file_path,parameters%sim_only)
        call GetString('input_folder',file_path,parameters%input_folder)
        call GetLogical('load_initial_guess',file_path,parameters%load_initial_guess)
        call GetInt('T_sim',file_path,parameters%T_sim,10)
        call GetInt('N_sim',file_path,parameters%N_sim,10)

        call GetReal('acc_int',file_path,parameters%acc_int,0.1_dp)
        call GetInt('acc_samples',file_path,parameters%acc_samples,10)
        call GetInt('acc_method',file_path,parameters%acc_method,1)

        call GetReal('perm_eps',file_path,parameters%perm_eps,0.01_dp)
        call GetReal('disc_share',file_path,parameters%disc_share,0.05_dp)
        call GetInt('T_sim_max',file_path,parameters%T_sim_max,10)
        call GetInt('RNG_seed',file_path,parameters%RNG_seed,0)
        call GetRealMat('mass',file_path,parameters%mass)
        call GetRealMat('l_max',file_path,parameters%l_max)
        call GetReal('beta',file_path,parameters%beta,0.95_dp)
        call GetRealMat('theta_0',file_path,parameters%theta_0)
        call GetRealMat('xi_theta',file_path,parameters%xi_theta)
        call GetReal('k_gov',file_path,parameters%k_gov,0.05_dp)
        call GetRealMat('S',file_path,parameters%S)
        call GetRealMat('P',file_path,parameters%P)
        call GetInt('s_init_index',file_path,parameters%s_init_index,1)
        call GetRealMat('b_min',file_path,parameters%b_min)
        call GetRealMat('b_max',file_path,parameters%b_max)
        call GetReal('B_min_kappa',file_path,parameters%B_min_kappa,0.1_dp)
        call GetReal('B_min_number',file_path,parameters%B_min_number,0.1_dp)
        call GetLogical('B_min_constant',file_path,parameters%B_min_constant)
        call GetRealMat('b_init',file_path,parameters%b_init)
        call GetRealMat('alpha',file_path,parameters%alpha)
        call GetReal('gamma_par',file_path,parameters%gamma_par,2.0_dp)
        call GetReal('A_par',file_path,parameters%A_par,1.0_dp)
        call GetReal('B_par',file_path,parameters%B_par,1.5_dp)
        call GetLogical('interior_solution_only',file_path,parameters%interior_solution_only)
        call GetLogical('one_fc_suff',file_path,parameters%one_fc_suff)
        call GetReal('stopping_rule_quantile',file_path,parameters%stopping_rule_quantile,1.0_dp)
        call GetReal('stopping_rule_threshold',file_path,parameters%stopping_rule_threshold,0.0000001_dp)
        call GetInt('VFI_max_iter',file_path,parameters%VFI_max_iter,1)
        call GetInt('VFI_glob_opt_iter',file_path,parameters%VFI_glob_opt_iter,0)
        call GetLogical('assume_interior_optimality',file_path,parameters%assume_interior_optimality)
        call GetInt('LFFC_max_ind',file_path,parameters%LFFC_max_ind,50)
        call GetReal('LFFC_min_share',file_path,parameters%LFFC_min_share,0.01_dp)
        call GetReal('LFFC_max_share',file_path,parameters%LFFC_max_share,0.9_dp)
        call GetInt('CTS_split_times',file_path,parameters%CTS_split_times,0)
        call GetInt('CTS_gridpoints_floor',file_path,parameters%CTS_gridpoints_floor,5)
        call GetInt('N_a',file_path,parameters%N_a,10)
        call GetInt('N_rho',file_path,parameters%N_rho,10)
        call GetInt('N_t',file_path,parameters%N_t,10)
        call GetReal('rho_min',file_path,parameters%rho_min,0.1_dp)
        call GetReal('rho_max',file_path,parameters%rho_max,2.0_dp)
        call GetRealMat('a_min',file_path,parameters%a_min)
        call GetRealMat('a_max',file_path,parameters%a_max)
        call GetInt('VFI_interpolation_mode',file_path,parameters%VFI_interpolation_mode,1)
        call GetReal('rho_search_mid',file_path,parameters%rho_search_mid,0.5_dp)
        call GetReal('rho_search_max_diff',file_path,parameters%rho_search_max_diff,0.5_dp)
        call GetReal('rho_search_floor',file_path,parameters%rho_search_floor,0.1_dp)
        call GetRealMat('a_search_min',file_path,parameters%a_search_min)
        call GetRealMat('a_search_max',file_path,parameters%a_search_max)
        call GetInt('N_search_a',file_path,parameters%N_search_a,10)
        call GetInt('N_search_rho',file_path,parameters%N_search_rho,10)
        call GetLogical('check_corners_only',file_path,parameters%check_corners_only)
        call GetLogical('stop_when_found',file_path,parameters%stop_when_found)
        call GetInt('save_iter_multiple',file_path,parameters%save_iter_multiple,1000)

        call GetInt('crit_How',file_path,parameters%crit_How,10)
        call GetRealMat('a_How',file_path,parameters%a_How)
        call GetIntMat('b_How',file_path,parameters%b_How)
        call GetReal('c_How',file_path,parameters%c_How,0.1_dp)

        call GetLogical('CTS_C_interp',file_path,parameters%CTS_C_interp)
        call GetLogical('save_C_guess',file_path,parameters%save_C_guess)

        call GetInt('acc_attmax1',file_path,parameters%acc_attmax1,100)
        call GetInt('acc_attmax2',file_path,parameters%acc_attmax2,10)
        call GetInt('acc_maxiter',file_path,parameters%acc_maxiter,500)
        call GetInt('LFFC_attmax',file_path,parameters%LFFC_attmax,10)
        call GetInt('LFFC_maxiter',file_path,parameters%LFFC_maxiter,200)
        call GetInt('VFI_opt_maxiter',file_path,parameters%VFI_opt_maxiter,200)
        call GetInt('sim_maxiter',file_path,parameters%sim_maxiter,200)
        call GetInt('sim_attmax1',file_path,parameters%sim_attmax1,100)
        call GetInt('sim_attmax2',file_path,parameters%sim_attmax2,10)
        call GetInt('opt_subroutine',file_path,parameters%opt_subroutine,2)

        call GetReal('rhoend',file_path,parameters%rhoend,0.000001_dp)
        call GetReal('e1',file_path,parameters%e1,0.000000001_dp)
        call GetReal('e2',file_path,parameters%e2,0.0000000001_dp)
        call GetReal('sim_rhoend',file_path,parameters%sim_rhoend,0.000001_dp)
        call GetReal('sim_e1',file_path,parameters%sim_e1,0.000000001_dp)
        call GetReal('sim_e2',file_path,parameters%sim_e2,0.0000000001_dp)


        do i=1,M_par !index of the row that we stack
            do j=1,M_par !stack it M times on top of each other
                parameters%P_onerow(i,j,:) = parameters%P(i,:)
            end do
        end do
        !relative growth rate
        parameters%rel_growth = (1.0_dp + parameters%xi_theta(1,2))/(1.0_dp + parameters%xi_theta(1,1))

end subroutine

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


        
        if(critical_error_found) then
                call double_output(folder_name,'Critical error found when checking parameters...')
                error stop
        end if
        

end subroutine check_parameters

end module mod_par
