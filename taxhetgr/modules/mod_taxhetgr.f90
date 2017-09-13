module mod_taxhetgr
!This module contains type definitions, subroutines, etc. which arewritten specifically for this project.
!(The division is not quite clear-cut of course, it is just a rough guideline)

use mod_types
use mod_par
use mod_parameters_hard
use mod_utility
use mod_tools

!USE nag_library
!Use nag_library, Only: e04jyf, nag_wp

implicit none

!Initial guess is composed of 2 types. The first type, c_guess, contains the guess at some binding constraint (or none).
!Type c_guess_all contains all the guesses at some gridpoint ('all' in the sense that these are guesses for
!all binding constraints or no constraint). It makes sense to split these because the type C_guess_all will
!be allocated memory depending on how much is needed (how many feasible guesses were find at that grid point).


!Type containing initial guess for consumption in solving the maximization problem at a given iteration.
!c_guess contains the actual guess (element (m,i) is consumption of agent i if the shock realization
!is m)., what_binds contains information on what constraint binds for this initial guess.
!what_binds = (i,j,k) means that agent i's constraint with index j bounds in state k, where j=-1 if lower bound
!on asset holdings binds, 1 if the upper bounds binds, and 0 if no constraint binds (interior solution)
type c_guess
    real(dp), dimension(M_par,I_par) :: c_guess
    integer, dimension(3) :: what_binds = [-1,-1,-1]
    !The default value is important in the program - where no feasible choice found,
    !what_binds is kept at this default value.
end type c_guess

!type c_guess_all contains, at a gridpoint, initial guesses for all binding constraints
!where an initial guess was found (including no constraint binding). Those constraints
!where nothing was found will be filled with special value (-1)
type c_guess_all
    type(c_guess), dimension(1+I_par*M_par) :: c_guess_all
end type c_guess_all

!type shared_data_c_guess contains two allocatable arrays, which can be used to share
!c_guess (typically one per gridpoint) between images. This has been separated because
!of issues with sharing coarrays of types which contain allocatable arrays of types -
!this issue is compiler-specific.

type shared_data_c_guess_type
    real(dp), allocatable, dimension(:,:,:) :: c !consumption allocation (last index notes grid point)
    integer, allocatable, dimension(:,:) :: wb !what binds (last index notes grid point)
end type shared_data_c_guess_type


!type grids_type contains grids. It is written for 2 agents only because
!the addition of third player would make the problem computationally infeasible
!using currently available hardware (going from 5 state variables to 7)
type grids_type
    real(dp),allocatable,dimension(:) :: A1_gr !asset holdings of agent 1
    real(dp),allocatable,dimension(:) :: A2_gr !asset holdings of agent 2
    real(dp),allocatable,dimension(:) :: rho_gr !u1/u2
    real(dp),allocatable,dimension(:) :: t_gr !relative productivity
end type grids_type

!type FP_result_type contains all sorts of results of the first-period problem
type FP_result_type
    real(dp), dimension(1,I_par) :: c,l,b_prime,a_prime,u_c,u_l
    real(dp) :: rho_prime !will need to increase dimension for I>2
    real(dp) :: t_prime !relative productivity of agent 2
    integer :: s_ind_init !initial index - this is usually drawn from a stationary distribution
    !(so that we don't get a solution dependent on initial state which can imply large quantitative
    !trend in government expenditure in expectation).
end type FP_result_type

!Type sim_result type contains
type sim_result_type
    !First index is time period, 2nd is agent
    real(dp), allocatable, dimension(:,:) :: c_sim,l_sim,b_sim,b_prime_sim,a_prime_sim,u_c_sim,u_l_sim,theta_lvl_sim
    !Note: there is no a_prime and rho_prime in the first period.
    real(dp), allocatable, dimension(:) :: rho_prime_sim
    real(dp), allocatable, dimension(:) :: t_prime_sim,t_sim,R_sim,tau_sim,g_lvl_sim
!    real(dp), allocatable, dimension(:) :: loss_sim !Extension.
    !Add stuff like Gini?
    integer, allocatable, dimension(:) :: s_ind_sim
    real(dp), allocatable, dimension(:) :: s_real_sim
    real(dp), allocatable, dimension(:) :: gini_income_sim
    real(dp), allocatable, dimension(:) :: gini_prod_sim
    real(dp), allocatable, dimension(:) :: gini_allincome_sim
    real(dp), allocatable, dimension(:) :: G_to_Y_sim
    real(dp), allocatable, dimension(:) :: Y_sim,B_gov_sim,B_gov_to_GDP_sim

    real(dp) :: a_slack_all_avg !average 'slackness' - basically low values (particularly negative) mean
    !that we have an issue in the sense that the a_prime state is close to bounds often, or even
    !gets outside of the bounds. The lower this is the worse the solution is expected to be and
    !this can then be used to eliminate 'outliers'.
end type sim_result_type

type sim_result_all_type
    type(sim_result_type), dimension(:), allocatable :: SRA !sim_result_all
    real(dp), dimension(:), allocatable :: tau_resid_last !residual in tax policy (relative to average) in the last period
    !(mreasure of variation between samples)
    real(dp), dimension(:), allocatable :: tau_resid_MAD !Max abs difference b/w residuals
    !(measure of variance withing a sample)
    real(dp), dimension(:), allocatable :: tau_resid_var !same as above but sample variance

    real(dp) :: avg_tau_MAD !max abs difference in tau in the sample average series (measure of
    !importance of trend but it might be better to compute this using a deterministic model
    !because the average will be also affected by uncertainty)

    real(dp) :: avg_tau_var !variance of average (across samples) tau.

    !Cross-sectional variance and MAD of residuals in the last period. This is intersting because
    !in a sense the residuals in the last period provide an upper bound on the (cumulative)
    !influence of shocks, and the cross-sectional variance of these thus provides
    !a certain measure of influence of shocks
    real(dp) :: tau_resid_CS_var,tau_resid_CS_MAD !CS - cross-sectional

    real(dp) :: tau_first,tau_last,tau_lastminfirst !average tax in period 1,N_sim, and period N_sim minus
    !period 1. - it's good summary of results to plot tau_first,tau_last as a function of alpha or something
    !in a single plot, using different colours!
    !And also plot tau_lastminfirst (this is essentially the same information but presented more clearly
    !to illustrate the difference in growth rate).

    !Same thing for transfers even though this is less interesting
    real(dp) :: trans_first,trans_last,trans_lastminfirst

    !relative significance of shocks:
    !1) within sample (relative to trend) - in terms of MAD and var ratios
        !avg(tau_resid_MAD)/avg_tau_MAD
        real(dp) :: stoch_metric_wi_MAD !stochasticity metric within sample (Max Abs Diff)
        !avg(tau_resid_var)/avg_tau_var
        real(dp) :: stoch_metric_wi_var !stochasticity metric within sample (variance)

    !2) across samples (cross-sectional variance of tau residuals in last period relative to trend
        !variation metric
        !tau_resid_CS_MAD/avg_tau_MAD
        real(dp) :: stoch_metric_bw_MAD !stochasticity metric between samples (Max Abs Diff)
        !tau_resid_CS_var/avg_tau_var
        real(dp) :: stoch_metric_bw_var !stochasticity metric between samples (variance)

    !average production for given shock realization - this is useful for calibration, such as when
    !we want to set the magnitude of shocks such that production drops by a certain percentage in
    !a recession (when g low)
    real(dp), dimension(M_par) :: avg_Y

    !Also save the 0.05 and 0.95 quantile of realizations of state variables. This is important for detecting
    !issues where the bounds of grid are wrongly chosen (the mean could still be within the bounds but
    !for some particular simulations (shock realizations the bounds could be hit or we could
    !be getting really close to them.
    !These definitions assume I=2, change later if needed!!!
    real(dp), dimension(:), allocatable :: rho_prime_05, rho_prime_95
    real(dp), dimension(:), allocatable :: a_prime_05, a_prime_95

end type sim_result_all_type


!Definition of operators for sim_result_type
interface operator(+)
    module procedure sim_plus_sim
end interface

interface operator(-)
    module procedure sim_min_sim
end interface

contains

!Definitions of functions used for overloading operators (such as for adding 2 sim_series, etc.)

!Adding simulation results (sr1, sr2)
!function sim_plus_sim(sr1,sr2)
function sim_plus_sim(sr1,sr2)
    type(sim_result_type) :: sim_plus_sim
    type(sim_result_type), intent(in) :: sr1,sr2

    type(par) :: pars_comm
    common /pars_comm/ pars_comm !common block parameters

    !Maybe there is a way to do this in a nicer way? This looks terrible and needs to be changed
    !if the definition of sim_result_type is changed...
    !If I used a subroutine and intent(inout), I wouldn't need to bother with the allocation,
    !but then I wouldn't be able to use operator overloading because for that I need a function...

    call sim_result_alloc(sim_plus_sim,pars_comm,1)

    sim_plus_sim%c_sim = sr1%c_sim + sr2%c_sim
    sim_plus_sim%l_sim = sr1%l_sim + sr2%l_sim
    sim_plus_sim%b_sim = sr1%b_sim + sr2%b_sim
    sim_plus_sim%b_prime_sim = sr1%b_prime_sim + sr2%b_prime_sim
!    sim_plus_sim%b_net_sim = sr1%b_net_sim + sr2%b_net_sim
!    sim_plus_sim%b_net_prime_sim = sr1%b_net_prime_sim + sr2%b_net_prime_sim
    sim_plus_sim%a_prime_sim = sr1%a_prime_sim + sr2%a_prime_sim
    sim_plus_sim%u_c_sim = sr1%u_c_sim + sr2%u_c_sim
    sim_plus_sim%u_l_sim = sr1%u_l_sim + sr2%u_l_sim
    sim_plus_sim%theta_lvl_sim = sr1%theta_lvl_sim + sr2%theta_lvl_sim
    sim_plus_sim%rho_prime_sim = sr1%rho_prime_sim + sr2%rho_prime_sim
    sim_plus_sim%t_prime_sim = sr1%t_prime_sim + sr2%t_prime_sim
    sim_plus_sim%t_sim = sr1%t_sim + sr2%t_sim
    sim_plus_sim%R_sim = sr1%R_sim + sr2%R_sim
    sim_plus_sim%tau_sim = sr1%tau_sim + sr2%tau_sim
    sim_plus_sim%g_lvl_sim = sr1%g_lvl_sim + sr2%g_lvl_sim
    sim_plus_sim%s_ind_sim = sr1%s_ind_sim + sr2%s_ind_sim
    sim_plus_sim%s_real_sim = sr1%s_real_sim + sr2%s_real_sim
    sim_plus_sim%gini_income_sim = sr1%gini_income_sim + sr2%gini_income_sim
    sim_plus_sim%gini_prod_sim = sr1%gini_prod_sim + sr2%gini_prod_sim
    sim_plus_sim%gini_allincome_sim = sr1%gini_allincome_sim + sr2%gini_allincome_sim
    sim_plus_sim%G_to_Y_sim = sr1%G_to_Y_sim + sr2%G_to_Y_sim
    sim_plus_sim%Y_sim = sr1%Y_sim + sr2%Y_sim
    sim_plus_sim%B_gov_sim = sr1%B_gov_sim + sr2%B_gov_sim
    sim_plus_sim%B_gov_to_GDP_sim = sr1%B_gov_to_GDP_sim + sr2%B_gov_to_GDP_sim
!    sim_plus_sim%trans_sim = sr1%trans_sim + sr2%trans_sim

end function sim_plus_sim

!Subtracting simulation results (sr1, sr2)
!function sim_min_sim(sr1,sr2)
function sim_min_sim(sr1,sr2)
    type(sim_result_type) :: sim_min_sim
    type(sim_result_type), intent(in) :: sr1,sr2

    type(par) :: pars_comm
    common /pars_comm/ pars_comm !common block parameters

    call sim_result_alloc(sim_min_sim,pars_comm,1)

    sim_min_sim%c_sim = sr1%c_sim - sr2%c_sim
    sim_min_sim%l_sim = sr1%l_sim - sr2%l_sim
    sim_min_sim%b_sim = sr1%b_sim - sr2%b_sim
    sim_min_sim%b_prime_sim = sr1%b_prime_sim - sr2%b_prime_sim
!    sim_min_sim%b_net_sim = sr1%b_net_sim - sr2%b_net_sim
!    sim_min_sim%b_net_prime_sim = sr1%b_net_prime_sim - sr2%b_net_prime_sim
    sim_min_sim%a_prime_sim = sr1%a_prime_sim - sr2%a_prime_sim
    sim_min_sim%u_c_sim = sr1%u_c_sim - sr2%u_c_sim
    sim_min_sim%u_l_sim = sr1%u_l_sim - sr2%u_l_sim
    sim_min_sim%theta_lvl_sim = sr1%theta_lvl_sim - sr2%theta_lvl_sim
    sim_min_sim%rho_prime_sim = sr1%rho_prime_sim - sr2%rho_prime_sim
    sim_min_sim%t_prime_sim = sr1%t_prime_sim - sr2%t_prime_sim
    sim_min_sim%t_sim = sr1%t_sim - sr2%t_sim
    sim_min_sim%R_sim = sr1%R_sim - sr2%R_sim
    sim_min_sim%tau_sim = sr1%tau_sim - sr2%tau_sim
    sim_min_sim%g_lvl_sim = sr1%g_lvl_sim - sr2%g_lvl_sim
    sim_min_sim%s_ind_sim = sr1%s_ind_sim - sr2%s_ind_sim
    sim_min_sim%s_real_sim = sr1%s_real_sim - sr2%s_real_sim
    sim_min_sim%gini_income_sim = sr1%gini_income_sim - sr2%gini_income_sim
    sim_min_sim%gini_prod_sim = sr1%gini_prod_sim - sr2%gini_prod_sim
    sim_min_sim%gini_allincome_sim = sr1%gini_allincome_sim - sr2%gini_allincome_sim
    sim_min_sim%G_to_Y_sim = sr1%G_to_Y_sim - sr2%G_to_Y_sim
    sim_min_sim%Y_sim = sr1%Y_sim - sr2%Y_sim
    sim_min_sim%B_gov_sim = sr1%B_gov_sim - sr2%B_gov_sim
    sim_min_sim%B_gov_to_GDP_sim = sr1%B_gov_to_GDP_sim - sr2%B_gov_to_GDP_sim
!    sim_min_sim%trans_sim = sr1%trans_sim - sr2%trans_sim

end function sim_min_sim

!Multiply simulation by (real) scalar
!function sim_multi(sr,multi)
function sim_multi(sr,multi)
    type(sim_result_type) :: sim_multi
    type(sim_result_type), intent(in) :: sr
    real(dp), intent(in) :: multi

    type(par) :: pars_comm
    common /pars_comm/ pars_comm !common block parameters

    call sim_result_alloc(sim_multi,pars_comm,1)

    sim_multi%c_sim = sr%c_sim * multi
    sim_multi%l_sim = sr%l_sim * multi
    sim_multi%b_sim = sr%b_sim * multi
    sim_multi%b_prime_sim = sr%b_prime_sim * multi
!    sim_multi%b_net_sim = sr%b_net_sim * multi
!    sim_multi%b_net_prime_sim = sr%b_net_prime_sim * multi
    sim_multi%a_prime_sim = sr%a_prime_sim * multi
    sim_multi%u_c_sim = sr%u_c_sim * multi
    sim_multi%u_l_sim = sr%u_l_sim * multi
    sim_multi%theta_lvl_sim = sr%theta_lvl_sim * multi
    sim_multi%rho_prime_sim = sr%rho_prime_sim * multi
    sim_multi%t_prime_sim = sr%t_prime_sim * multi
    sim_multi%t_sim = sr%t_sim * multi
    sim_multi%R_sim = sr%R_sim * multi
    sim_multi%tau_sim = sr%tau_sim * multi
    sim_multi%g_lvl_sim = sr%g_lvl_sim * multi
    sim_multi%s_ind_sim = sr%s_ind_sim * multi
    sim_multi%s_real_sim = sr%s_real_sim * multi
    sim_multi%gini_income_sim = sr%gini_income_sim * multi
    sim_multi%gini_prod_sim = sr%gini_prod_sim * multi
    sim_multi%gini_allincome_sim = sr%gini_allincome_sim * multi
    sim_multi%G_to_Y_sim = sr%G_to_Y_sim * multi
    sim_multi%Y_sim = sr%Y_sim * multi
    sim_multi%B_gov_sim = sr%B_gov_sim * multi
    sim_multi%B_gov_to_GDP_sim = sr%B_gov_to_GDP_sim * multi
!    sim_multi%trans_sim = sr%trans_sim * multi
end function sim_multi


!Divide simulation by (real) scalar
!function sim_div(sr,div)
function sim_div(sr,div)
    type(sim_result_type) :: sim_div
    type(sim_result_type), intent(in) :: sr
    real(dp), intent(in) :: div

    real(dp) :: div_inverse

    type(par) :: pars_comm
    common /pars_comm/ pars_comm !common block parameters

    call sim_result_alloc(sim_div,pars_comm,1)

    div_inverse = 1/div

    sim_div = sim_multi(sr,div_inverse)

end function sim_div

!Subroutine b_init_perm generates a random permutation of the initial conditions.
!This subroutine does not check that the permutation doesn't violate borrowing constraints,
!so care must be taken in this case (eps shouldn't be too large and pars%b_init
!should not be too close to constraints).
subroutine b_init_permutate(b_init,b_init_perm,perm_eps)
    real(dp), dimension(1,I_par), intent(in) :: b_init !initial conditions as given in parameter file
    real(dp), dimension(1,I_par), intent(out) :: b_init_perm !permutated initial conditions
    real(dp) :: perm_eps

    real(dp) :: rnd_01 !random number drawn from 0,1 interval (uniform distribution)

    call random_number(rnd_01)
    b_init_perm(1,1) = b_init(1,1) + perm_eps*(2.0_dp*rnd_01 - 1.0_dp)

    call random_number(rnd_01)
    b_init_perm(1,2) = b_init(1,2) + perm_eps*(2.0_dp*rnd_01 - 1.0_dp)


end subroutine b_init_permutate

!Function b_resid computes residual = b_prime - b_constr.
function b_resid(c_im,b_constr,i_binds,m_binds,pars,cons_mat,g_lvl,theta_lvl,a_lvl,P_onerow)
real(dp) :: b_resid
real(dp), intent(in) :: c_im,b_constr
integer, intent(in) :: i_binds,m_binds
type(par), intent(in) :: pars
real(dp), dimension(M_par,I_par), intent(in) :: cons_mat !incomplete consumption matrix
real(dp), dimension(M_par,1), intent(in) :: g_lvl !government expenditure in levels (not normalized by productivity)
real(dp), dimension(1,I_par), intent(in) :: theta_lvl
real(dp), dimension(I_par), intent(in) :: a_lvl
real(dp), dimension(M_par,M_par), intent(in) :: P_onerow

real(dp), dimension(M_par,I_par) :: cons_mat_comp !complete cons. mat
real(dp), dimension(M_par,I_par) :: l !labour
real(dp), dimension(M_par,I_par) :: a_prime,b_prime !next-period asset holdings

!complete the consumption matrix using the guessed element
cons_mat_comp = cons_mat
cons_mat_comp(M_par,I_par) = c_im
!!!Make sure to use the complete consumption matrix in the remainder of this function

!Get labour supply
call lab(cons_mat_comp,l,theta_lvl,pars%mass,g_lvl,pars%Gamma_par)

!Get next-period assets
call ass(a_lvl,cons_mat_comp,l,pars%A_par,pars%B_par,pars%Gamma_par,pars%beta,P_onerow,a_prime)

!Get asset holdings
call get_b_prime(b_prime,a_prime,cons_mat_comp,pars)

!Compute the residual. Also normalize by productivity - by current period
!productivity since that's what the constraint is in terms of!!!
!This is different to previous versions of the paper, leads to
!more efficiency (because we don't need to compute next-period productivity) to
!look for feasible bounds, at least as long as we don't use state variables
!normalized by next-period productivity (which we don't in early versions
!of this program).
b_resid = b_prime(i_binds,m_binds)/theta_lvl(1,i_binds) - b_constr

end function b_resid


!Simple subroutine which converts a matrix of guesses and what_binds point to
!type c_guess.
subroutine c_guess_gen(c_guess_mat,what_binds,c_guess_type)
    real(sp), dimension(M_par,I_par), intent(in) :: c_guess_mat
    integer, dimension(3), intent(in) :: what_binds
    type(c_guess), intent(out) :: c_guess_type

    c_guess_type%c_guess = c_guess_mat
    c_guess_type%what_binds = what_binds

end subroutine c_guess_gen

!Subroutine c_inc_mtv transforms the incomplete consumption matrix (with element
!M_par,I_par missing) into a row form, so it can be used by some NAG functions
!for maximization. It first uses reshape on the first M_par - 1 rows to
!put them into the first (M_par-1) * I_par elements of the vector. Then it
!takes the last row of the matrix (M_par) and puts its first I_par - 1
!at the end of the vector.
subroutine c_inc_mtv(c_inc_mat,c_inc_vec)
    real(dp),dimension(M_par,I_par), intent(in) :: c_inc_mat
    real(dp), dimension(M_par*I_par - 1), intent(out) :: c_inc_vec

    c_inc_vec(1:(M_par-1)*I_par) = reshape(c_inc_mat(1:M_par-1,1:I_par),[(M_par-1)*I_par])
    c_inc_vec((M_par-1)*I_par+1:ncvar_par) = c_inc_mat(M_par,1:I_par-1)

    !ncvar_par = M_par*I_par - 1 is the number of choice variables
end subroutine c_inc_mtv

!subroutine c_inc_vtm is pretty much just an inverse of subroutine c_inc_mtv
subroutine c_inc_vtm(c_inc_vec,c_inc_mat)
    real(dp), dimension(M_par*I_par - 1), intent(in) :: c_inc_vec
    real(dp),dimension(M_par,I_par), intent(inout) :: c_inc_mat !inout because
    !the last element of the matrix is undefined so we want to leave it at the previous
    !value

    c_inc_mat(1:M_par-1,1:I_par) = reshape(c_inc_vec(1:(M_par-1)*I_par),[M_par-1,I_par])
    c_inc_mat(M_par,1:I_par-1) = c_inc_vec((M_par-1)*I_par+1:ncvar_par)

end subroutine c_inc_vtm

!In the latest version of the program, interior solutions only are used (but the program still
!contains some references to borrowing constraints which were in an earlier version. We need
!to make sure that the computations are run with interior_solution_only = .true.).

!Subroutine check_constr checks whether constraints are violated. It also computes a convex
!cost of violating constraints (so far this is equal to 0)
subroutine check_constr(c,l,b_prime,b_min_mat,b_max_mat,l_max_mat,theta_lvl,a_prime,a_min_mat,a_max_mat&
,what_binds,rho_min_vec,rho_max_vec,constr_ok,loss,A_par,rofv,neg_loss)
    real(dp), dimension(M_par,I_par), intent(in) :: c,l,b_prime,a_prime
    real(dp), dimension(M_par,I_par), intent(in) :: l_max_mat
    real(dp), dimension(M_par,I_par), intent(in) :: b_min_mat, b_max_mat, a_min_mat,a_max_mat
    integer, dimension(3), intent(in) :: what_binds
    real(dp), dimension(M_par,1), intent(in) :: rho_min_vec,rho_max_vec
    real(dp), dimension(1,I_par), intent(in) :: theta_lvl
    real(dp), intent(in) :: A_par
    logical, intent(in) :: neg_loss !If this is true then loss is negative for points which strictly satisfy
    !all constraints, and it is 'more negative' the further from bounds the point is. This is used to find
    !the 'most feasible' point in the sense of laying in the middle of some large convex feasible set (we need
    !to do this because of nonconvexity of feasible sets).


    logical, intent(out) :: constr_ok !this will be true if all constraints are satisfied
    real(dp),intent(out) :: loss !loss stemming from violating constraints.

    logical, intent(in) :: rofv !return on first violation
    !If this is true, once a constraint is violated, the control is returned, and
    !we don't finish checking the other constraints. This should be false when looking for
    !feasible choices until we find one feasible choice, then it should be false
    !(no point looking for not strictly feasible choices if we already have a strictly
    !feasible one). In the maximization (computing continuation return), this should
    !be false.
    !This speeds up computation substantially.

    real(dp), dimension(M_par,I_par) :: b_prime_norm !normalized asset holdings (divided by next-period productivity)

    real(dp), dimension(M_par,1) :: ones = 1.0_dp

    real(dp), dimension(M_par,I_par) :: u_c !marginal utility of consumption
    real(dp), dimension(M_par,1) :: rho_prime

    real(dp), dimension(M_par,I_par) :: b_min_mat_local, b_max_mat_local !Local copies: These will
    !be used to check constraints (will be essentially the same as the matrices passed to this subroutine
    !possibly with one constraint replaced by a large constant).

    real(dp) :: loss1a,loss2a,loss3a !levels of loss from least severe problem to the most severe

!Get a copy of the constraints
b_min_mat_local = b_min_mat
b_max_mat_local = b_max_mat

!The function checks all constraints - even if one is violated we still go through till the end
!and add a cost of violation which should make it unappealing for the program to choose such point.
!
!There will be two loss conditions: a strict one for violating a fundamental constraint, and
!a milder one for violating constrain on marginal-utility adjusted

!parameters of loss function directly here
!later move them to parameter file.
!Functional form of loss: quadratic for every violation of a constraint. This can be experimented on.
loss1a = 100.0_dp
loss2a = 1000.0_dp
loss3a = 100000.0_dp

!
!Also the functional forms could change, not just the parameters.

!initialize constraints_ok and loss
constr_ok = .true.
loss = 0.0_dp

!On violating a critical constraint, return control to the calling process once the first violation
!of a constraint is detected. With less critical constraints introduce loss later (use merge),
!also see notes on my tablet.

!(1) Check whether that consumption is positive
if(any(c<=0)) then
    constr_ok = .false.
    loss = 1000000000_dp !Negative consumption is a huge problem, therefore a large constant penalty is applied
    !so it will never be chosen in the optimization
    !even at the risk of introducing a discontinuity. In the program we should never actually get here as the maximizations
    !should be constrained and happening over a positive region of consumption.
    return
    !return always if this is violated, negative consumption never allowed. Perhaps experiment
    !with rofv and continuous loss here if the discontinuity causes problems.
 end if

!(2) Check that labour lies between 0 and l_max. If Labour is negative, apply
!large constant loss as in case of consumption. If it exceeds l_max, apply a convex cost.
if(any(l<0)) then
    constr_ok = .false.
    if(rofv) return !return before even computing loss because if rofv is true,
    !we only care about strict feasibility and not about value of loss function.

    !sum of squared violations multiplied by a parameter.
    loss = loss + loss3a*sum(min(l,0.0_dp)**2)

end if

if(any(l>l_max_mat)) then
    constr_ok = .false.
    if(rofv) return !return before even computing loss because if rofv is true,
    !we only care about strict feasibility and not about value of loss function.

!sum of squared violations multiplied by a parameter.
loss = loss + loss2a*sum(max(l-l_max_mat,0.0_dp)**2)

end if

!(3) Check whether borrowing constraints are not violated. This condition is perhaps a bit less
!important than the conditions on consumption and labour so we don't need to be quite as strict

!theta_prime_mat = matmul(ones,theta_prime)
!write(*,*) theta_prime_mat
b_prime_norm = b_prime/matmul(ones,theta_lvl) !Normalize by current productivity since that's what
!the constraint for b is in terms of (in recent version of the paper)


!If none of the constraints binds (what_binds(1) = 0), then check all constraints.
!Otherwise do not check the one that is supposed to bind. This is distinguished at the outset
!to increase efficiency.

!(In the current version of TaxHetGr interior solutions only should be used)

if(what_binds(1)==0) then !Interior solution, check all constraints.

    if(any(b_prime_norm<b_min_mat_local)) then
        constr_ok = .false.
        if(rofv) return
        loss = loss + loss2a*sum(max(b_min_mat_local-b_prime_norm,0.0_dp)**2)
    end if
else !Not an interior solution. In this case we ignore the one constraint which is supposed to bind.
    !Replace the value of that constrain in the constraint matrix to arbitrarily large or low value
    !so that it never binds. Just make sure we don't use the overwritten constraints anywehre else.
    if(what_binds(2)==-1) then !lower constraint of some agent binds
        b_min_mat_local(what_binds(3),what_binds(1)) = -1000000000000
            !Remember that we replaced the constraint - when computing the loss that one number in the matrix
            !needs to be set to 0.
        if(any(b_prime_norm>b_max_mat_local)) then
            constr_ok = .false.
            if(rofv) return
            loss = loss + loss2a*sum(max(b_prime_norm-b_max_mat_local,0.0_dp)**2)
        end if
        if(any(b_prime_norm<b_min_mat_local)) then
            constr_ok = .false.
            if(rofv) return
            loss = loss + loss2a*sum(max(b_min_mat_local-b_prime_norm,0.0_dp)**2)
        end if

    end if
    if(what_binds(2)==1) then !upper constraint of some agent binds
        b_max_mat_local(what_binds(3),what_binds(1)) = 1000000000000

        if(any(b_prime_norm>b_max_mat_local)) then
            constr_ok = .false.
            if(rofv) return
            loss = loss + loss2a*sum(max(b_prime_norm-b_max_mat_local,0.0_dp)**2)
        end if
        if(any(b_prime_norm<b_min_mat_local)) then
            constr_ok = .false.
            if(rofv) return
            loss = loss + loss2a*sum(max(b_min_mat_local-b_prime_norm,0.0_dp)**2)
        end if
    end if
!This looks messy but efficiency is important here, really can't afford any do cycles so deep in the program.

end if

!(4) Checking constraints on a_prime.

!Here loss function will be particularly important (no jumps)
if(any(a_prime>a_max_mat)) then
    constr_ok = .false.
    if(rofv) return
    loss = loss + loss1a*sum(max(a_prime-a_max_mat,0.0_dp)**2)
end if
if(any(a_prime<a_min_mat)) then
    constr_ok = .false.
    if(rofv) return
    loss = loss + loss1a*sum(max(a_min_mat-a_prime,0.0_dp)**2)
end if

!(EFE): This subroutine should already take rho as input -> sometimes we compute it twice (first when
!checking cosntraints and also when computing the next-period return).
!
!Maybe this is also the case of another of the variables passed here (such as b_?)

!(5) Check constraints on rho
!First need to compute the marginal utility of consumption from c and l
call util_c(c,u_c,A_par)

rho_prime = reshape(u_c(:,1)/u_c(:,2),[M_par,1])

if(any(rho_prime>rho_max_vec)) then
    constr_ok = .false.
    if(rofv) return
    loss = loss + loss1a*sum(max(rho_prime-rho_max_vec,0.0_dp)**2)
end if
if(any(rho_prime<rho_min_vec)) then
    constr_ok = .false.
    if(rofv) return
    loss = loss + loss1a*sum(max(rho_min_vec-rho_prime,0.0_dp)**2)
end if

!If neg_loss (negative loss) is true, we subtract a number for points which lie strictly inside the feasible set.
!This of course cannot be used in value iteration or simulation because it would distort the returns. However, this
!is useful if we look for a feasible choice and we want to find something strictly inside a feasible set (because
!the loss function otherwise stops decreasing once we reach a boundary of the feasible set and there is no
!benefit in reaching a strictly feasible point).
!
!The measure used here is euclidean metric, i.e., distance from the relevant borders (in one direction only).
!This way the return from interior feasibility is approximately linear so it should not be beneficial
!to get a point strictly feasible in some dimensions while infeasible in others, because the
!penalizations there are much higher.
if(neg_loss) then
    !first add distance from zero for consumption (other things being equal, we want to favour points with high consumption)
    loss = loss - sqrt(sum(max(c,0.0_dp)**2))

    !Do the same adjustment for lower and upper bounds for labour
    loss = loss - sqrt(sum(max(l,0.0_dp)**2))

    !And for state variables (this has no economic meaning, just that we want to start looking for
    !the solution as far from bounds as possible. If it's suboptimal it will drift from this region anyway).
    !But if we start looking for a solution in a 'corner' region of the state space, we might get
    !anomalous behaviour.

    !For a_prime scale the adjustment by the range. Otherwise the terms here will grow linearly with the width
    !of feasible set for MU adjusted asset holdings. It could then possibly outweigh loss and make us choose
    !points which are in the middle of feasible set for a_prime but are infeasible with respect to
    !other constraints. Also it would favour feasibility in sense of a_prime relative to feasibility in sense of c or l.
    !The normalization assumes that all elements of a_max_mat - a_min_mat are positive. If not we would
    !get all sorts of problems, crashes, and non-sensical results elsewhere.
    loss = loss - sqrt(sum(max((a_prime-a_min_mat)/(a_max_mat-a_min_mat),0.0_dp)**2))
    loss = loss - sqrt(sum(min((a_prime-a_max_mat)/(a_max_mat-a_min_mat),0.0_dp)**2))

    !And lastly do the same for rho_prime. This is relevant only in case of stochastic calibration or
    !corner solutions. In the case of deterministic calibration and interior solutions the value
    !of rho is fixed so it will be just a constant term in objective function of a LFFC loss minimization
    !subroutine.
    loss = loss - sqrt(sum(max((rho_prime-rho_min_vec)/(rho_max_vec-rho_min_vec),0.0_dp)**2))
    loss = loss - sqrt(sum(min((rho_prime-rho_max_vec)/(rho_max_vec-rho_min_vec),0.0_dp)**2))

end if

end subroutine check_constr

!Subroutine CPC_ret computes (-1)*( current + continuation return), given current state
!and the incomplete consumption matrix (in vector form).
!
!It is written in such form that it can be used by NAG Fortran Subroutine E04JYF
!and possibly other subroutines. For this reason, many of the input parameters
!are passed to this function using common block.
!

subroutine CPC_ret(n,xc,fc,iuser,ruser)
        real (dp), intent (out) :: fc
        integer, intent (in) :: n
        real (dp), intent (inout) :: ruser(*)
        real (dp), intent (in) :: xc(n)
        integer, intent (inout) :: iuser(*) !In this instance, it is used to pass what_binds

        real(dp), dimension(:,:,:,:,:), pointer :: V_old_pntr !pointer to V_old, passed using common block

        !The following variables are passed to the function using a common block
        real(dp), dimension(:), pointer :: a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr
        real(dp), dimension(:,:), pointer :: S_pntr,P_pntr, mass_pntr
        !and transition matrix.
        real(dp), pointer :: A_par_pntr,B_par_pntr,Gamma_par_pntr
        type(par), pointer :: pars_pntr

        !Values of state variables as saved in the grid (hence the suffix st).
        real(dp), dimension(I_par) :: a_st
        real(dp), dimension(I_par-1) :: rho_st
        real(dp) :: t_st,g_st
        integer :: g_st_ind
        real(dp), dimension(1,I_par) :: theta_lvl

        real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat
        real(dp), dimension(M_par,1) :: rho_min_vec,rho_max_vec

        real(dp), dimension(M_par,I_par) :: g_lvl_mat
        real(dp), dimension(M_par,M_par) :: P_onerow

        !Local variables
        real(dp), dimension(M_par,I_par) :: c_inc_mat = 0.0_dp !incomplete consumption matrix
        !Labour
        real(dp), dimension(M_par,I_par) :: lab_mat = 0.0_dp
        !next-period states
        real(dp) :: t_st_pr
        real(dp), dimension(M_par,I_par) :: a_prime,b_prime

        logical :: c_im_found = .false., constr_ok

        integer :: g_ind !index of shock realization

        real(dp) :: loss !loss from violating constraints.

        real(dp), dimension(M_par,1) :: conval_all !continuation value in all states

        !rho_prime and u_c
        real(dp), dimension(M_par,I_par) :: u_c !marginal utility of consumption
        real(dp), dimension(M_par,1) :: rho_prime


        common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
        !to some subroutines.
        common /pars_pointer/ pars_pntr

        common /cpc_ret_bloc/ a_st,theta_lvl,rho_st,t_st,g_st,g_st_ind
        common /cpc_ret_bloc2/ g_lvl_mat, P_onerow,l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec

        common /comm1/ a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr,A_par_pntr,B_par_pntr,&
        Gamma_par_pntr,S_pntr,P_pntr,mass_pntr

        common /commloss/ loss

        !Get the incomplete consumption matrix from the input vector (save it as c_inc_mat)
        call c_inc_vtm(xc,c_inc_mat)


        !Depending on what_binds (stored in iuser), compute the missing element of the consumption
        !matrix. If no element is found, which should not happen, use 0.0001_dp, which should
        !result in sufficiently low utility so that this point will never be optimal. The
        !discontinuity stemming from this step could have severe implications for
        !finding optimum though... Maybe think of something else. Such as choosing
        !last-period choice as optimum, etc.
        !OMG -> here I need the whole pars object... This will be problematic for generalizations
        !in the future when it may contain allocatable objects. Will need to
        !rewrite get_c_im and all subroutines it uses s.t. here we only need the relevant
        !bits of pars.
        call get_c_im(c_inc_mat,iuser,a_st,rho_st,theta_lvl,g_lvl_mat(:,1),g_st_ind,pars_pntr,P_onerow,c_im_found)
        !(from this point on, c_inc_mat is actually the complete consumption plan)

        if(.not. c_im_found) then
            fc = 10000.0_dp !significant discontinuity, hopefully will almost never happen
            return
        end if



        !Compute Labour (for that we need the state passed in a common block)
        call lab(c_inc_mat,lab_mat,theta_lvl,mass_pntr,g_lvl_mat(:,1),Gamma_par_pntr)

        !Now that we know consumption and labour plan, we can compute the current (expected) return
        !Note: will need to pass mass separately and get rid of pars_pntr.
        fc = SWF(c_inc_mat,lab_mat,g_st_ind,P_pntr,mass_pntr,pars_pntr%alpha,A_par_pntr,B_par_pntr,Gamma_par_pntr)


        !Compute next-period states (need to compute a_prime and theta_prime)
        t_st_pr = t_prime(t_st,t_gr_pntr(pars_pntr%N_t))


        !(EFE) : Get rid of this and use rel_prod_prime instead, more efficient here...

        call ass(a_st,c_inc_mat,lab_mat,pars_pntr%A_par,pars_pntr%B_par,pars_pntr%Gamma_par,&
        pars_pntr%beta,P_onerow,a_prime)

        call get_b_prime(b_prime,a_prime,c_inc_mat,pars_pntr)

        constr_ok = .false.
        !Check constraints (Get a loss value). All of the matrices/vectors for checking constraints need to
        !be precomputed and passed to this thing by common block.

!        !debug
!        if(this_image()==1) then
!        write(*,*) 'bounds within CPC_ret 1st'
!        write(*,*) 'b_min_mat = ',b_min_mat
!        write(*,*) 'b_max_mat = ',b_max_mat
!        write(*,*) 'l_max_mat = ',l_max_mat
!        write(*,*) 'a_min_mat = ',a_min_mat
!        write(*,*) 'a_max_mat = ',a_max_mat
!        write(*,*) '_______________________'
!        end if

        call check_constr(c_inc_mat,lab_mat,b_prime,b_min_mat,b_max_mat,l_max_mat,theta_lvl,a_prime,a_min_mat&
        ,a_max_mat,iuser,rho_min_vec,rho_max_vec,constr_ok,loss,pars_pntr%A_par,.false.,.false.)
        !Now add discounted expected continuation value: This will be computed by function dECV

!        !debug
!        if(this_image()==1) then
!        write(*,*) 'bounds within CPC_ret 2nd'
!        write(*,*) 'b_min_mat = ',b_min_mat
!        write(*,*) 'b_max_mat = ',b_max_mat
!        write(*,*) 'l_max_mat = ',l_max_mat
!        write(*,*) 'a_min_mat = ',a_min_mat
!        write(*,*) 'a_max_mat = ',a_max_mat
!        write(*,*) '_______________________'
!        end if

        !Get rho_prime. This is not very efficient because function check_constr computes it too.
        !(EFE) In the future create a subroutine for computing rho_prime, then pass it to
        !check_constr already as an argument. This is quite important because this subroutine
        !is used extremely often so efficiency is important.
        call util_c(c_inc_mat,u_c,pars_pntr%A_par)
        rho_prime = reshape(u_c(:,1)/u_c(:,2),[M_par,1])

        !If the choice is not strictly feasible, force it to be within bounds (extrapolation
        !is very bad, usually leads to divergence) if there is no feasible choice.
        !There is potential gain for efficiency here if we check feasibility constraint by
        !constraint and force only the problematic point to be in the bounds.
        !Only need to do this for a_prime and rho_prime, theta and shock realization will
        !be inside grid by construction!
        if(.not. constr_ok) then
            a_prime=max(a_min_mat,a_prime)
            a_prime=min(a_max_mat,a_prime)
            rho_prime=max(rho_min_vec,rho_prime)
            rho_prime=min(rho_max_vec,rho_prime)
        end if


        !Once we have the point forced within the grid boundaries, get the continuation
        !value by interpolation. We need to get a separate continuation value for all
        !states.
        do g_ind = 1,M_par
            call interp_V(V_old_pntr,[a_prime(g_ind,1),a_prime(g_ind,2),rho_prime(g_ind,1),t_st_pr]&
            ,g_ind,conval_all(g_ind,1),a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr,size(a1_gr_pntr),&
            size(rho_gr_pntr),size(t_gr_pntr))
        end do

!        write(*,*) 'Debug in subroutine CPC_ret:'
!        write(*,*) 'cons = ',c_inc_mat
!        write(*,*) 'current return = ',fc
!        write(*,*) 'conval = ',conval_all
!        write(*,*) 'loss = ',loss
!        write(*,*) '______________________'

        !Now we have continuation value for every possible shock realization. Add the discounted
        !expectation of this term to the current return. (g_st_ind is the last shock realization)

        fc = fc + pars_pntr%beta*maxval(matmul(P_pntr(g_st_ind,:),conval_all))
        !maxval is there simply to convert (1,1) array to scalar. Maybe this can be done in a better way...

        !Subtract loss stemming from violating some constraints.
        fc = (fc - loss)*(-1.0_dp)
        !We need negative return here!!!



end subroutine CPC_ret

!Subroutine CPC_ret2 is a copy of CPC_ret, the only difference being the presence of one
!more output. This is defined for use in a different NAg subroutine than CPC_ret.
subroutine CPC_ret2(n,xc,fc,iuser,ruser,inform)
        real (dp), intent (out) :: fc
        integer, intent (in) :: n
        real (dp), intent (inout) :: ruser(*)
        real (dp), intent (in) :: xc(n)
        integer, intent (inout) :: iuser(*) !In this instance, it is used to pass what_binds
        integer, intent(out) :: inform

        real(dp), dimension(:,:,:,:,:), pointer :: V_old_pntr !pointer to V_old, passed using common block

        !The following variables are passed to the function using a common block
        real(dp), dimension(:), pointer :: a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr
        real(dp), dimension(:,:), pointer :: S_pntr,P_pntr, mass_pntr
        !and transition matrix.
        real(dp), pointer :: A_par_pntr,B_par_pntr,Gamma_par_pntr
        type(par), pointer :: pars_pntr

        !Values of state variables as saved in the grid (hence the suffix st).
        real(dp), dimension(I_par) :: a_st
        real(dp), dimension(I_par-1) :: rho_st
        real(dp) :: t_st,g_st
        integer :: g_st_ind
        real(dp), dimension(1,I_par) :: theta_lvl

        real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat
        real(dp), dimension(M_par,1) :: rho_min_vec,rho_max_vec

        real(dp), dimension(M_par,I_par) :: g_lvl_mat
        real(dp), dimension(M_par,M_par) :: P_onerow

        !Local variables
        real(dp), dimension(M_par,I_par) :: c_inc_mat = 0.0_dp !incomplete consumption matrix
        !Labour
        real(dp), dimension(M_par,I_par) :: lab_mat = 0.0_dp
        !next-period states
        real(dp) :: t_st_pr
        real(dp), dimension(M_par,I_par) :: a_prime,b_prime

        logical :: c_im_found = .false., constr_ok

        integer :: g_ind !index of shock realization

        real(dp) :: loss !loss from violating constraints.

        real(dp), dimension(M_par,1) :: conval_all !continuation value in all states

        !rho_prime and u_c
        real(dp), dimension(M_par,I_par) :: u_c !marginal utility of consumption
        real(dp), dimension(M_par,1) :: rho_prime


        common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
        !to some subroutines.
        common /pars_pointer/ pars_pntr

        common /cpc_ret_bloc/ a_st,theta_lvl,rho_st,t_st,g_st,g_st_ind
        common /cpc_ret_bloc2/ g_lvl_mat, P_onerow,l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec

        common /comm1/ a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr,A_par_pntr,B_par_pntr,&
        Gamma_par_pntr,S_pntr,P_pntr,mass_pntr

        common /commloss/ loss

        !Negative value lead to termination of maximization procedure.
        inform = 0

        !Get the incomplete consumption matrix from the input vector (save it as c_inc_mat)
        call c_inc_vtm(xc,c_inc_mat)


        call get_c_im(c_inc_mat,iuser,a_st,rho_st,theta_lvl,g_lvl_mat(:,1),g_st_ind,pars_pntr,P_onerow,c_im_found)
        !(from this point on, c_inc_mat is actually the complete consumption plan)
        if(.not. c_im_found) then
            fc = 10000.0_dp !significant discontinuity, hopefully will almost never happen
        end if

        !Compute Labour (for that we need the state passed in a common block)
        call lab(c_inc_mat,lab_mat,theta_lvl,mass_pntr,g_lvl_mat(:,1),Gamma_par_pntr)

        !Now that we know consumption and labour plan, we can compute the current (expected) return
        !Note: will need to pass mass separately and get rid of pars_pntr.
        fc = SWF(c_inc_mat,lab_mat,g_st_ind,P_pntr,mass_pntr,pars_pntr%alpha,A_par_pntr,B_par_pntr,Gamma_par_pntr)

        !Compute next-period states
        t_st_pr = t_prime(t_st,t_gr_pntr(pars_pntr%N_t))

        call ass(a_st,c_inc_mat,lab_mat,pars_pntr%A_par,pars_pntr%B_par,pars_pntr%Gamma_par,&
        pars_pntr%beta,P_onerow,a_prime)

        call get_b_prime(b_prime,a_prime,c_inc_mat,pars_pntr)

        constr_ok = .false.
        !Check constraints (Get a loss value). All of the matrices/vectors for checking constraints need to
        !be precomputed and passed to this thing by common block.
        call check_constr(c_inc_mat,lab_mat,b_prime,b_min_mat,b_max_mat,l_max_mat,theta_lvl,a_prime,a_min_mat&
        ,a_max_mat,iuser,rho_min_vec,rho_max_vec,constr_ok,loss,pars_pntr%A_par,.false.,.false.)
        !Now add discounted expected continuation value: This will be computed by function dECV

        !Get rho_prime. This is not very efficient because function check_constr computes it too.
        !(EFE) In the future create a subroutine for computing rho_prime, then pass it to
        !check_constr already as an argument. This is quite important because this subroutine
        !is used extremely often so efficiency is important.
        call util_c(c_inc_mat,u_c,pars_pntr%A_par)
        rho_prime = reshape(u_c(:,1)/u_c(:,2),[M_par,1])

        !If the choice is not strictly feasible, force it to be within bounds (extrapolation
        !is very bad, usually leads to divergence) if there is no feasible choice.
        !There is potential gain for efficiency here if we check feasibility constraint by
        !constraint and force only the problematic point to be in the bounds.
        !Only need to do this for a_prime and rho_prime, theta and shock realization will
        !be inside grid by construction!
        if(.not. constr_ok) then
            a_prime=max(a_min_mat,a_prime)
            a_prime=min(a_max_mat,a_prime)
            rho_prime=max(rho_min_vec,rho_prime)
            rho_prime=min(rho_max_vec,rho_prime)
        end if


        !Once we have the point forced within the grid boundaries, get the continuation
        !value by interpolation. We need to get a separate continuation value for all
        !states.
        do g_ind = 1,M_par
            call interp_V(V_old_pntr,[a_prime(g_ind,1),a_prime(g_ind,2),rho_prime(g_ind,1),t_st_pr]&
            ,g_ind,conval_all(g_ind,1),a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr,size(a1_gr_pntr),&
            size(rho_gr_pntr),size(t_gr_pntr))
        end do

        !Now we have continuation value for every possible shock realization. Add the discounted
        !expectation of this term to the current return. (g_st_ind is the last shock realization)

        fc = fc + pars_pntr%beta*maxval(matmul(P_pntr(g_st_ind,:),conval_all))
        !maxval is there simply to convert (1,1) array to scalar. Maybe this can be done in a better way...

        !Subtract loss stemming from violating some constraints.
        fc = (fc - loss)*(-1.0_dp)
        !We need negative return here!!!


end subroutine CPC_ret2

!Subroutine CPC_ret3 is a copy of CPC_ret, the only difference being the presence 2 more inputs
!This is defined for use in a different NAg subroutine than CPC_ret.

subroutine CPC_ret3(n,xc,fc,nstate,iuser,ruser,inform)
        real (dp), intent (out) :: fc
        integer, intent (in) :: n
        integer, intent(in) :: nstate !=1 if the subroutine is called for the first time,
        !not actually used anywhere atm but required by NAg.
        real (dp), intent (inout) :: ruser(*)
        real (dp), intent (in) :: xc(n)
        integer, intent (inout) :: iuser(*) !In this instance, it is used to pass what_binds
        integer, intent(out) :: inform

        real(dp), dimension(:,:,:,:,:), pointer :: V_old_pntr !pointer to V_old, passed using common block

        !The following variables are passed to the function using a common block
        real(dp), dimension(:), pointer :: a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr
        real(dp), dimension(:,:), pointer :: S_pntr,P_pntr, mass_pntr
        !and transition matrix.
        real(dp), pointer :: A_par_pntr,B_par_pntr,Gamma_par_pntr
        type(par), pointer :: pars_pntr

        !Values of state variables as saved in the grid (hence the suffix st).
        real(dp), dimension(I_par) :: a_st
        real(dp), dimension(I_par-1) :: rho_st
        real(dp) :: t_st,g_st
        integer :: g_st_ind
        real(dp), dimension(1,I_par) :: theta_lvl

        real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat
        real(dp), dimension(M_par,1) :: rho_min_vec,rho_max_vec

        real(dp), dimension(M_par,I_par) :: g_lvl_mat
        real(dp), dimension(M_par,M_par) :: P_onerow

        !Local variables
        real(dp), dimension(M_par,I_par) :: c_inc_mat = 0.0_dp !incomplete consumption matrix
        !Labour
        real(dp), dimension(M_par,I_par) :: lab_mat = 0.0_dp
        !next-period states
        real(dp) :: t_st_pr
        real(dp), dimension(M_par,I_par) :: a_prime,b_prime

        logical :: c_im_found = .false., constr_ok

        integer :: g_ind !index of shock realization

        real(dp) :: loss !loss from violating constraints.

        real(dp), dimension(M_par,1) :: conval_all !continuation value in all states

        !rho_prime and u_c
        real(dp), dimension(M_par,I_par) :: u_c !marginal utility of consumption
        real(dp), dimension(M_par,1) :: rho_prime


        common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
        !to some subroutines.
        common /pars_pointer/ pars_pntr

        common /cpc_ret_bloc/ a_st,theta_lvl,rho_st,t_st,g_st,g_st_ind
        common /cpc_ret_bloc2/ g_lvl_mat, P_onerow,l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec

        common /comm1/ a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr,A_par_pntr,B_par_pntr,&
        Gamma_par_pntr,S_pntr,P_pntr,mass_pntr

        common /commloss/ loss

        !Negative value lead to termination of maximization procedure.
        inform = 0

        !Get the incomplete consumption matrix from the input vector (save it as c_inc_mat)
        call c_inc_vtm(xc,c_inc_mat)


        call get_c_im(c_inc_mat,iuser,a_st,rho_st,theta_lvl,g_lvl_mat(:,1),g_st_ind,pars_pntr,P_onerow,c_im_found)
        !(from this point on, c_inc_mat is actually the complete consumption plan)
        if(.not. c_im_found) then
            fc = 10000.0_dp !significant discontinuity, hopefully will almost never happen
        end if

        !Compute Labour (for that we need the state passed in a common block)
        call lab(c_inc_mat,lab_mat,theta_lvl,mass_pntr,g_lvl_mat(:,1),Gamma_par_pntr)

        !Now that we know consumption and labour plan, we can compute the current (expected) return
        !Note: will need to pass mass separately and get rid of pars_pntr.
        fc = SWF(c_inc_mat,lab_mat,g_st_ind,P_pntr,mass_pntr,pars_pntr%alpha,A_par_pntr,B_par_pntr,Gamma_par_pntr)

        !Compute next-period states
        t_st_pr = t_prime(t_st,t_gr_pntr(pars_pntr%N_t))

        call ass(a_st,c_inc_mat,lab_mat,pars_pntr%A_par,pars_pntr%B_par,pars_pntr%Gamma_par,&
        pars_pntr%beta,P_onerow,a_prime)

        call get_b_prime(b_prime,a_prime,c_inc_mat,pars_pntr)

        constr_ok = .false.
        !Check constraints (Get a loss value). All of the matrices/vectors for checking constraints need to
        !be precomputed and passed to this thing by common block.
        call check_constr(c_inc_mat,lab_mat,b_prime,b_min_mat,b_max_mat,l_max_mat,theta_lvl,a_prime,a_min_mat&
        ,a_max_mat,iuser,rho_min_vec,rho_max_vec,constr_ok,loss,pars_pntr%A_par,.false.,.false.)
        !Now add discounted expected continuation value: This will be computed by function dECV

        !Get rho_prime. This is not very efficient because function check_constr computes it too.
        !(EFE) In the future create a subroutine for computing rho_prime, then pass it to
        !check_constr already as an argument. This is quite important because this subroutine
        !is used extremely often so efficiency is important.
        call util_c(c_inc_mat,u_c,pars_pntr%A_par)
        rho_prime = reshape(u_c(:,1)/u_c(:,2),[M_par,1])

        !If the choice is not strictly feasible, force it to be within bounds (extrapolation
        !is very bad, usually leads to divergence) if there is no feasible choice.
        !There is potential gain for efficiency here if we check feasibility constraint by
        !constraint and force only the problematic point to be in the bounds.
        !Only need to do this for a_prime and rho_prime, theta and shock realization will
        !be inside grid by construction!
        if(.not. constr_ok) then
            a_prime=max(a_min_mat,a_prime)
            a_prime=min(a_max_mat,a_prime)
            rho_prime=max(rho_min_vec,rho_prime)
            rho_prime=min(rho_max_vec,rho_prime)
        end if


        !Once we have the point forced within the grid boundaries, get the continuation
        !value by interpolation. We need to get a separate continuation value for all
        !states.
        do g_ind = 1,M_par
            call interp_V(V_old_pntr,[a_prime(g_ind,1),a_prime(g_ind,2),rho_prime(g_ind,1),t_st_pr]&
            ,g_ind,conval_all(g_ind,1),a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr,size(a1_gr_pntr),&
            size(rho_gr_pntr),size(t_gr_pntr))
        end do

        !Now we have continuation value for every possible shock realization. Add the discounted
        !expectation of this term to the current return. (g_st_ind is the last shock realization)

        fc = fc + pars_pntr%beta*maxval(matmul(P_pntr(g_st_ind,:),conval_all))
        !maxval is there simply to convert (1,1) array to scalar. Maybe this can be done in a better way...

        !Subtract loss stemming from violating some constraints.
        fc = (fc - loss)*(-1.0_dp)


end subroutine CPC_ret3

!The following is a one-dimensional version of CPC which is used when M=1 and I=2.
!(some multidimensional optimization subroutines do not work properly in this case)
subroutine CPC_ret_1d(xc,fc,iuser,ruser)
        real (dp), intent (out) :: fc
        real (dp), intent (inout) :: ruser(*)
        real (dp), intent (in) :: xc
        integer, intent (inout) :: iuser(*) !In this instance, it is used to pass what_binds

        real(dp), dimension(:,:,:,:,:), pointer :: V_old_pntr !pointer to V_old, passed using common block

        !The following variables are passed to the function using a common block
        real(dp), dimension(:), pointer :: a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr
        real(dp), dimension(:,:), pointer :: S_pntr,P_pntr, mass_pntr
        !and transition matrix.
        real(dp), pointer :: A_par_pntr,B_par_pntr,Gamma_par_pntr
        type(par), pointer :: pars_pntr

        !Values of state variables as saved in the grid (hence the suffix st).
        real(dp), dimension(I_par) :: a_st
        real(dp), dimension(I_par-1) :: rho_st
        real(dp) :: t_st,g_st
        integer :: g_st_ind
        real(dp), dimension(1,I_par) :: theta_lvl

        real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat
        real(dp), dimension(M_par,1) :: rho_min_vec,rho_max_vec

        real(dp), dimension(M_par,I_par) :: g_lvl_mat
        real(dp), dimension(M_par,M_par) :: P_onerow

        !Local variables
        real(dp), dimension(M_par,I_par) :: c_inc_mat = 0.0_dp !incomplete consumption matrix
        !Labour
        real(dp), dimension(M_par,I_par) :: lab_mat = 0.0_dp
        !next-period states
        real(dp) :: t_st_pr
        real(dp), dimension(M_par,I_par) :: a_prime,b_prime

        logical :: c_im_found = .false., constr_ok

        integer :: g_ind !index of shock realization

        real(dp) :: loss !loss from violating constraints.

        real(dp), dimension(M_par,1) :: conval_all !continuation value in all states

        !rho_prime and u_c
        real(dp), dimension(M_par,I_par) :: u_c !marginal utility of consumption
        real(dp), dimension(M_par,1) :: rho_prime


        common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
        !to some subroutines.
        common /pars_pointer/ pars_pntr

        common /cpc_ret_bloc/ a_st,theta_lvl,rho_st,t_st,g_st,g_st_ind
        common /cpc_ret_bloc2/ g_lvl_mat, P_onerow,l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec

        common /comm1/ a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr,A_par_pntr,B_par_pntr,&
        Gamma_par_pntr,S_pntr,P_pntr,mass_pntr

        common /commloss/ loss

        !Get the incomplete consumption matrix from the input vector (save it as c_inc_mat)
        !call c_inc_vtm(xc,c_inc_mat)
        c_inc_mat = xc

        call get_c_im(c_inc_mat,iuser,a_st,rho_st,theta_lvl,g_lvl_mat(:,1),g_st_ind,pars_pntr,P_onerow,c_im_found)
        !(from this point on, c_inc_mat is actually the complete consumption plan)
        if(.not. c_im_found) then
            fc = 10000.0_dp !significant discontinuity, hopefully will almost never happen
            return
        end if

        !Compute Labour (for that we need the state passed in a common block)
        call lab(c_inc_mat,lab_mat,theta_lvl,mass_pntr,g_lvl_mat(:,1),Gamma_par_pntr)

        !Now that we know consumption and labour plan, we can compute the current (expected) return
        !Note: will need to pass mass separately and get rid of pars_pntr.
        fc = SWF(c_inc_mat,lab_mat,g_st_ind,P_pntr,mass_pntr,pars_pntr%alpha,A_par_pntr,B_par_pntr,Gamma_par_pntr)

        !Compute next-period states
        t_st_pr = t_prime(t_st,t_gr_pntr(pars_pntr%N_t))

        call ass(a_st,c_inc_mat,lab_mat,pars_pntr%A_par,pars_pntr%B_par,pars_pntr%Gamma_par,&
        pars_pntr%beta,P_onerow,a_prime)

        call get_b_prime(b_prime,a_prime,c_inc_mat,pars_pntr)

        constr_ok = .false.
        !Check constraints (Get a loss value). All of the matrices/vectors for checking constraints need to
        !be precomputed and passed to this thing by common block.
        call check_constr(c_inc_mat,lab_mat,b_prime,b_min_mat,b_max_mat,l_max_mat,theta_lvl,a_prime,a_min_mat&
        ,a_max_mat,iuser,rho_min_vec,rho_max_vec,constr_ok,loss,pars_pntr%A_par,.false.,.false.)
        !Now add discounted expected continuation value: This will be computed by function dECV

        !Get rho_prime. This is not very efficient because function check_constr computes it too.
        !(EFE) In the future create a subroutine for computing rho_prime, then pass it to
        !check_constr already as an argument. This is quite important because this subroutine
        !is used extremely often so efficiency is important.
        call util_c(c_inc_mat,u_c,pars_pntr%A_par)
        rho_prime = reshape(u_c(:,1)/u_c(:,2),[M_par,1])


        if(.not. constr_ok) then
            a_prime=max(a_min_mat,a_prime)
            a_prime=min(a_max_mat,a_prime)
            rho_prime=max(rho_min_vec,rho_prime)
            rho_prime=min(rho_max_vec,rho_prime)
        end if


        !Once we have the point forced within the grid boundaries, get the continuation
        !value by interpolation. We need to get a separate continuation value for all
        !states.
        do g_ind = 1,M_par
            call interp_V(V_old_pntr,[a_prime(g_ind,1),a_prime(g_ind,2),rho_prime(g_ind,1),t_st_pr]&
            ,g_ind,conval_all(g_ind,1),a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr,size(a1_gr_pntr),&
            size(rho_gr_pntr),size(t_gr_pntr))
        end do

        !Now we have continuation value for every possible shock realization. Add the discounted
        !expectation of this term to the current return. (g_st_ind is the last shock realization)

        fc = fc + pars_pntr%beta*maxval(matmul(P_pntr(g_st_ind,:),conval_all))
        !maxval is there simply to convert (1,1) array to scalar. Maybe this can be done in a better way...

        !Subtract loss stemming from violating some constraints.
        fc = (fc - loss)*(-1.0_dp)
        !We need negative return here!!!

end subroutine CPC_ret_1d


!Subroutine CTS_num_points computes the number of grid points at an iteration of CTS algorithm.
subroutine CTS_num_points(CTS_ind,pars,N_a,N_rho,N_t)
    integer, intent(in) :: CTS_ind
    type(par), intent(in) :: pars
    integer, intent(out) :: N_a,N_rho,N_t

    integer :: split = 0 !times that the number of gridpoints is split at this iteration of CTS algorithm

    if(pars%CTS_split_times>0) then
        split = pars%CTS_split_times + 1 - CTS_ind

        N_a = floor(pars%N_a * 0.5_dp**split);
        N_rho = floor(pars%N_rho * 0.5_dp**split);
        N_t = floor(pars%N_t * 0.5_dp**split);

        !If some resulting number of grid points is lower than the minimum number of grid points
        !per that variable, adjust this. However, do this only if pars%CTS_split_times > 0. If this is
        !zero it means that we are not using CTS algorithm and it makes no sense to use the floor.
        N_a=maxval([pars%CTS_gridpoints_floor,N_a])
        N_rho=maxval([pars%CTS_gridpoints_floor,N_rho])
        N_t=maxval([pars%CTS_gridpoints_floor,N_t])
    else !CTS algorithm not used, use the number of gridpoints given in parameter file
        N_a = pars%N_a
        N_rho = pars%N_rho
        N_t = pars%N_t
    end if

    !If static_environment is used, it means that there is no growth. In this case
    !set N_theta = 1.
    if(pars%static_environment) N_t = 1

end subroutine CTS_num_points


!subroutine gen_equi_grid generates an equispaced grid.
subroutine gen_equi_grid(grid,point1,pointN)
    real(dp), dimension(:) :: grid
    real(dp), intent(in) :: point1,pointN

    integer :: N
    real(dp) :: dist !distance b/w end and final point
    integer :: i

    !If point1 > pointN, it's a mistake (could allow this and reverse the
    !order of points but this will help reveal bugs)
    if(point1>=pointN) then
        write(*,*) 'Error: point1 >= pointN in subroutine gen_equi_grid.'
        error stop
    end if

    N = size(grid)

    dist = (pointN - point1)/real(N-1)

    grid(1) = point1

    do i=2,N-1
        grid(i) = grid(i-1) + dist
    end do
    grid(N) = pointN !do this separately to avoid problems caused by imprecision.

end subroutine gen_equi_grid


!Subroutine gen_WB_poss generates all possible combinations of what_binds.
subroutine gen_WB_poss(WB_poss)
    integer, dimension(1+I_par*M_par,3), intent(inout) :: WB_poss

    integer :: i,m, poss_ind !index for agent and shock realization, and possibility
    !The first element says whose constraint binds (if 0 no constraint binds),
    !the second says whether it is lower or upper bound (-1 for lower, 1 for upper),
    !the last says in what state the constraint binds. We stack
    !all possibilities into matrix where every row is a possible what_binds

    !First the case where no constraint binds
    WB_poss(1,:) = [0,0,0]
    poss_ind=2
    do i=1,I_par !cycle over agents
        do m=1,M_par !cycle over states
            !lower bound
            WB_poss(poss_ind,:) = [i,-1,m]
            poss_ind = poss_ind + 1
        end do
    end do



end subroutine gen_WB_poss

!This subroutine will compute b_prime in all states.
subroutine get_b_prime(b_prime,a_prime,c,pars)
   real(dp), dimension(M_par,I_par), intent(inout) :: b_prime
   real(dp), dimension(M_par,I_par), intent(in) :: a_prime
   real(dp), dimension(M_par,I_par), intent(in) :: c
   type(par), intent(in) :: pars

   real(dp), dimension(M_par,I_par) :: u_c
   !get marginal utility of consumption
   call util_c(c,u_c,pars%A_par)

    !remember: a_i =(u_i*b_i)/beta
   b_prime = pars%beta*a_prime/u_c

end subroutine get_b_prime


!subroutine get_c_guess_img computes an initial guess of consumption for all agents in all states
!at every gridpoint on which an image works. It also computes the share of nodes at which this is feasible.
!This result is saved in coarray share_feas_img.
!The feasible guesses are returned in vector of type c_guess_all (see its definition in this module).
subroutine get_c_guess_img(c_guess_img,min_ind_img,max_ind_img,ind_unfold_all,pars,grids,share_feas_img,&
avg_loss_img,LFFC_choices,num_fc,CTS_ind,C_pol_all,grids_tmp)
    real(dp), codimension[*] :: share_feas_img !The weighted share of points where a feasible c was found
    real(dp) :: avg_loss_img
    integer, intent(in) :: min_ind_img,max_ind_img
    type(par), codimension[*] :: pars !parameters. Not sure why this is a coarray here... Could be a coarray
    !only where we need to share data, here it could be just a standard type. Low priority issue though.
    type(c_guess_all), dimension(min_ind_img:max_ind_img), intent(out) :: c_guess_img
    integer, dimension(5,min_ind_img:max_ind_img) :: ind_unfold_all
    real(dp), dimension(:,:,:), intent(in) :: LFFC_choices
    integer, intent(in) :: num_fc
    type(grids_type), intent(in) :: grids
    integer, intent(in) :: CTS_ind
    real(dp), dimension(:,:,:,:,:,:,:), allocatable, intent(in) :: C_pol_all
    type(grids_type), intent(in) :: grids_tmp

    integer :: gp_ind !index for cycling over gridpoints. The main index from which we will get the indices
    !for the state variables from ind_unfold_all
    !Some of the following variables will become row vectors.
    !integer :: a1_ind,a2_ind,rho_ind,theta_ind,m_ind !- indices not needed to be kept separately?

    !Values of state variables as saved in the grid (hence the suffix st).
    real(dp), dimension(I_par) :: a_st
    real(dp), dimension(I_par-1) :: rho_st
    real(dp) :: t_st,g_st
    integer :: g_st_ind

    !Some variables used in calling subroutine LFFC
    real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat
    real(dp), dimension(M_par,1) :: rho_min_vec,rho_max_vec
    type(c_guess_all) :: c_guess_all_gp !output: all feasible choices at a grid point
    logical :: fc_found !true if a feasible choice found at this gridpoint
    !possibilities of what_binds
    integer, dimension(1+I_par*M_par,3) :: WB_poss

    !real(dp) :: rel_prod_min !minimum relative productivity. This will be the lowest relative productivity
    !in the grid for this variable. This is used in computing next-period states by LFFC.

    real(dp) :: loss_gp, loss_img !loss at gridpoint, loss at image (average weighted by number of
    !gridpoints at this image).


    if(I_par>2) then
        write(*,*) 'Error: subroutine get_c_guess_img assumes I_par = 2.'
        error stop
    end if


    !Get bounds in matrix/vector form to save some computation time later
    call prep_bounds(pars%l_max,pars%b_min,pars%b_max,pars%a_min,pars%a_max,pars%rho_min,pars%rho_max,&
    l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec)
    !Note: here we use bounds from parameters. In the algorithm searching for feasible bounds
    !we will use bounds depending on stage of the algorithm (bounds for a and rho will vary)

    !Generate possible what_binds
    call gen_WB_poss(WB_poss)


    !initialize average loss and feasible share
    avg_loss_img = 0.0_dp
    share_feas_img = 0.0_dp
    do gp_ind=min_ind_img,max_ind_img
        !find values of state variables and call subroutine LFFC which finds a feasible choice
        !given current state. Here just recover the raw states (as captured in the grid) and pass
        !them to LFFC.
        a_st(1) = grids%a1_gr(ind_unfold_all(1,gp_ind))
        a_st(2) = grids%a2_gr(ind_unfold_all(2,gp_ind))
        rho_st = grids%rho_gr(ind_unfold_all(3,gp_ind))
        t_st = grids%t_gr(ind_unfold_all(4,gp_ind))
        g_st = pars%S(1,ind_unfold_all(5,gp_ind))
        g_st_ind = ind_unfold_all(5,gp_ind)


        call LFFC(a_st,rho_st,t_st,g_st,g_st_ind,pars,b_min_mat,b_max_mat,l_max_mat,a_min_mat,a_max_mat,&
        rho_min_vec,rho_max_vec,LFFC_choices,num_fc,WB_poss,c_guess_all_gp,fc_found,loss_gp,CTS_ind,C_pol_all,grids_tmp)

        !save the guess
        c_guess_img(gp_ind) = c_guess_all_gp

        !update statistics on feasibility of shares, etc.
        avg_loss_img = avg_loss_img + loss_gp/real(max_ind_img-min_ind_img+1)

        if(fc_found) then
            share_feas_img = share_feas_img + 1.0_dp/real(max_ind_img-min_ind_img+1)
        end if

    end do

    !debug:
    !write(*,*) 'img ',this_image(),'share_feas_img = ',share_feas_img
end subroutine get_c_guess_img



!subroutine get_c_im computes element I,M of the consumption matrix, depending on what_binds,
!i.e., on which borrowing constraint binds (if any).
subroutine get_c_im(cons_mat,what_binds,a_st,rho_st,theta_lvl,g_lvl_vec,g_ind_st,pars,P_onerow,c_im_found)
    real(dp), dimension(M_par,I_par), intent(inout) :: cons_mat
    integer, dimension(3), intent(in) :: what_binds
    real(dp), dimension(I_par), intent(in) :: a_st
    real(dp), dimension(I_par-1), intent(in) :: rho_st
    real(dp), dimension(1,I_par), intent(in) :: theta_lvl
    real(dp), dimension(M_par,1), intent(in) :: g_lvl_vec !gov expenditure in levels in all states
    integer, intent(in) :: g_ind_st !index
    type(par), intent(in) :: pars
    real(dp), dimension(M_par,M_par), intent(in) :: P_onerow
    logical, intent(inout) :: c_im_found !when it enters it is false, change if something found

    integer :: constr_ag,constr_st !constrained player and state in which he is constrained.
    real(dp) :: constr_value !value of the constraint

    real(dp), dimension(M_par,I_par) :: b_prime,l

    real(dp) :: c_im_tmp

    real(dp) :: x_l,x_m,x_h,fx_l,fx_h,fx_m !used in bisection algorithm
    integer :: bisec_ind
    real(dp) :: bisec_tol

    bisec_tol = 0.000001_dp

    !If interior solution, use function CM.
    if(what_binds(1) == 0) then
        c_im_tmp = c2M_int(cons_mat,rho_st(1),pars%P,g_ind_st)
        cons_mat(M_par,I_par) = c_im_tmp
        if(c_im_tmp>0.0_dp) then
            c_im_found = .true.
        else
            !In previous versions this was absent so if (wrongly) c_im_found was true
            !on entry and no feasible choice was found, then the error was not detected,
            !and it looked as if we found non-negative c_IM.
            !Also we used to copy the valye found by cons_mat only if c_im_found was true.
            !Now we copy it anyway - it should not be used anywhere, but it is useful for debugging.
            c_im_found = .false.
        end if
        return
    end if



    if(what_binds(2) == -1) then !lower constraint
        constr_value = pars%b_min(1,what_binds(1))
    else !Just assume that what_binds was generated properly (saves time here)
        constr_value = pars%b_max(1,what_binds(1))
    end if


    !Bisection: (the following could be moved into a subroutine)

    !1) Initialize two points for c_{I,M}
    x_l = 0.001_dp
    x_h = (sum(pars%mass*theta_lvl*pars%l_max) - g_lvl_vec(what_binds(3),1))/pars%mass(1,what_binds(1));

    if(x_h<x_l) x_h=1.0_dp

    !First make sure that we have different signs at end-points. If not, return a non-sensical value
    !which will be rejected by the algorithm at later stages.
    fx_l = b_resid(x_l,constr_value,what_binds(1),what_binds(3),pars,cons_mat,g_lvl_vec,theta_lvl,a_st,P_onerow)
    fx_h = b_resid(x_h,constr_value,what_binds(1),what_binds(3),pars,cons_mat,g_lvl_vec,theta_lvl,a_st,P_onerow)

    !If not strictly less than 0, exit. Strictly speaking, equality
    !would mean that one of the endpoints exactly satisfies the equation but that will
    !never happen in practice
    if(fx_l*fx_h >= 0) then
        !exiting now implies that c_im_found will be .false.
        return
    end if

    !Now continue with the bisection algorithm.
    do bisec_ind = 1,100 !max number of iterations set here
        x_m = (x_l + x_h)/2.0_dp
        fx_m = b_resid(x_m,constr_value,what_binds(1),what_binds(3),pars,cons_mat,g_lvl_vec,theta_lvl,a_st,P_onerow)

        if(fx_m*fx_l < 0) then
            x_h = x_m
            fx_h = fx_m
        else
            x_l = x_m
            fx_l = fx_m
        end if

        if(abs(x_h-x_l) < bisec_tol) exit
    end do

    cons_mat(M_par,I_par) = (x_l + x_h)/2.0_dp
    c_im_found = .true.

end subroutine get_c_im

!subroutine for generating grids. 4 equispaced grids are generated for the first four states
!(a_1,a_2,rho,t).
!The grid for shock space is loaded from parameters and we keep it separate (this grid is
!never resized).
subroutine grids_gen(A1_gr,A2_gr,rho_gr,t_gr,N_a,N_rho,N_t,pars)
    real(dp), dimension(:), intent(out) :: A1_gr,A2_gr,rho_gr,t_gr
    integer, intent(in) :: N_a,N_rho,N_t !can't just load these from parameter type
    !because in CTS algorithm they vary.
    type(par) :: pars


    call gen_equi_grid(A1_gr,pars%a_min(1,1),pars%a_max(1,1))

    call gen_equi_grid(A2_gr,pars%a_min(1,2),pars%a_max(1,2))

    call gen_equi_grid(rho_gr,pars%rho_min,pars%rho_max)

    !Grid for 'time'
    if(N_T==1) then !this should happen only if static_environment == 1.
        T_gr(1) = 0.0_dp !Fix 'time' at t=0.
        return
    end if

    !if not N_t == 1, then the first element will be t=1, and the last one
    !t=T_sim_max. We never need the continuation value for t=0.
    call gen_equi_grid(T_gr,1.0_dp,real(pars%T_sim_max,dp))

end subroutine

!Subroutine interp_V returns interpolated value of value function at a point
!anywhere in the range of the grids.
!(It assumes monotonically increasing grids, as the subroutine locate_gr
!which is used makes that assumption).
!
!V is the value function (usually V_old), x is the point in state space
!at which we want to get the value function (a1,a2,rho,theta) - we never
!interpolate across shock space (s_ind is the fixed index(!) of shock realization),
!s_ind is index of shock realization, y is the interpolated value.
!The rest of the inputs are grids for the state variables, and number
!of elements of these grids (so we don't have to use function size() which
!may be slightly slower).
!
!The function also works for extrapolation. If N_theta is 1 (static environment),
!we use a separate formula for increased efficiency.
!
!The interpolation used here is trilinear or quadrilinear depending on whether
!theta_gr has more than one element (static or dynamic environment).
!
!The function also allows extrapolation and does not perform any checks
!as to whether we are on the grid (to increase efficiency).
subroutine interp_V(V,x,s_ind,y,a1_gr,a2_gr,rho_gr,t_gr,N_a,N_rho,N_t)
    real(dp), dimension(:,:,:,:,:), intent(in) :: V
    real(dp), dimension(4), intent(in) :: x
    integer, intent(in) :: s_ind
    real(dp), intent(out) :: y
    real(dp), dimension(:) ,intent(in) :: a1_gr,a2_gr,rho_gr,t_gr

    integer, intent(in) :: N_a,N_rho,N_t
    !Same number of gridpoints for agent 1 and agent 2 is assumed - this can be
    !generalized easily, but its difficult to see the benefit of doing it.
    integer :: j_a1,j_a2,j_rho,j_t !used for storing which points are to be used in interpolation.
    real(dp) :: t_a1,t_a2,t_rho,t_t ! weights in linear combinations t*x(j) + (1-t)*x(j+1)


    !This follows quite closely algorithm suggested in the Numerical recipes in Fortran book.
    !1) Find the points to be used in interpolation (for all variables, point with index j
    !is the point s.t. x lies between x(j) and x(j+1).
    call locate_gr(a1_gr,N_a,x(1),j_a1) !j_a1 will be the index of the grid point in a1_gr
    call locate_gr(a2_gr,N_a,x(2),j_a2)
    call locate_gr(rho_gr,N_rho,x(3),j_rho)


    !coefficients in linear combinations (value 0 is associated with the lower point)
    t_a1 = (a1_gr(j_a1+1)-x(1))/(a1_gr(j_a1+1)-a1_gr(j_a1))
    t_a2 = (a2_gr(j_a2+1)-x(2))/(a2_gr(j_a2+1)-a2_gr(j_a2))
    t_rho = (rho_gr(j_rho+1)-x(3))/(rho_gr(j_rho+1)-rho_gr(j_rho))

    !split The algorithm into 2 branches depending on whether we are also
    !interpolating over theta (static_environment = .false.) or not!
    if(N_t>1) then
        !4-dimensional interpolation
        call locate_gr(t_gr,N_t,x(4),j_t)
        t_t = (t_gr(j_t+1)-x(4))/(t_gr(j_t+1)-t_gr(j_t))

        !The following looks horrible but it should be ok (generated using a spreadsheet script)
        !(this subroutine is used potentially millions of times per one iteration in VFI algorithm
        !so maximum efficiency is important - no saving of intermediate results or calling
        !subroutines). Also, elements of V are accessed in the order in which they are stored
        !in memory - it could be the case that the compiler does optimization of this sort automatically
        !but it doesn't hurt to write it explicitly.
        y = V(j_a1,j_a2,j_rho,j_t,s_ind)*t_a1*t_a2*t_rho*t_t+&
            V(j_a1+1,j_a2,j_rho,j_t,s_ind)*(1-t_a1)*t_a2*t_rho*t_t+&
            V(j_a1,j_a2+1,j_rho,j_t,s_ind)*t_a1*(1-t_a2)*t_rho*t_t+&
            V(j_a1+1,j_a2+1,j_rho,j_t,s_ind)*(1-t_a1)*(1-t_a2)*t_rho*t_t+&
            V(j_a1,j_a2,j_rho+1,j_t,s_ind)*t_a1*t_a2*(1-t_rho)*t_t+&
            V(j_a1+1,j_a2,j_rho+1,j_t,s_ind)*(1-t_a1)*t_a2*(1-t_rho)*t_t+&
            V(j_a1,j_a2+1,j_rho+1,j_t,s_ind)*t_a1*(1-t_a2)*(1-t_rho)*t_t+&
            V(j_a1+1,j_a2+1,j_rho+1,j_t,s_ind)*(1-t_a1)*(1-t_a2)*(1-t_rho)*t_t+&
            V(j_a1,j_a2,j_rho,j_t+1,s_ind)*t_a1*t_a2*t_rho*(1-t_t)+&
            V(j_a1+1,j_a2,j_rho,j_t+1,s_ind)*(1-t_a1)*t_a2*t_rho*(1-t_t)+&
            V(j_a1,j_a2+1,j_rho,j_t+1,s_ind)*t_a1*(1-t_a2)*t_rho*(1-t_t)+&
            V(j_a1+1,j_a2+1,j_rho,j_t+1,s_ind)*(1-t_a1)*(1-t_a2)*t_rho*(1-t_t)+&
            V(j_a1,j_a2,j_rho+1,j_t+1,s_ind)*t_a1*t_a2*(1-t_rho)*(1-t_t)+&
            V(j_a1+1,j_a2,j_rho+1,j_t+1,s_ind)*(1-t_a1)*t_a2*(1-t_rho)*(1-t_t)+&
            V(j_a1,j_a2+1,j_rho+1,j_t+1,s_ind)*t_a1*(1-t_a2)*(1-t_rho)*(1-t_t)+&
            V(j_a1+1,j_a2+1,j_rho+1,j_t+1,s_ind)*(1-t_a1)*(1-t_a2)*(1-t_rho)*(1-t_t)


    else
        !3-dimensional interpolation (trilinear). In this case, the fourth index
        !is fixed, so we work with V(:,:,:,1,s_ind).
        y = V(j_a1,j_a2,j_rho,1,s_ind)*t_a1*t_a2*t_rho + &
            V(j_a1+1,j_a2,j_rho,1,s_ind)*(1-t_a1)*t_a2*t_rho +&
            V(j_a1,j_a2+1,j_rho,1,s_ind)*t_a1*(1-t_a2)*t_rho +&
            V(j_a1+1,j_a2+1,j_rho,1,s_ind)*(1-t_a1)*(1-t_a2)*t_rho +&
            V(j_a1,j_a2,j_rho+1,1,s_ind)*t_a1*t_a2*(1-t_rho) +&
            V(j_a1+1,j_a2,j_rho+1,1,s_ind)*(1-t_a1)*t_a2*(1-t_rho) +&
            V(j_a1,j_a2+1,j_rho+1,1,s_ind)*t_a1*(1-t_a2)*(1-t_rho) +&
            V(j_a1+1,j_a2+1,j_rho+1,1,s_ind)*(1-t_a1)*(1-t_a2)*(1-t_rho)

    end if

end subroutine interp_V

!subroutine interp_C serves for interpolating policy function - it is essentially pretty much the same
!as interp_V, the difference being that we get an interpolated consumption matrix (where interpolation
!is done for each element separately).
subroutine interp_C(C_pol_all,x,s_ind,y,a1_gr,a2_gr,rho_gr,t_gr,N_a,N_rho,N_t)
    real(dp), dimension(:,:,:,:,:,:,:), intent(in) :: C_pol_all !policy function
    real(dp), dimension(4), intent(in) :: x
    integer, intent(in) :: s_ind
    real(dp), intent(out), dimension(M_par,I_par) :: y !interpolated choice matrix (state-contingent cons plan forall agents)
    real(dp), dimension(:), intent(in) :: a1_gr,a2_gr,rho_gr,t_gr

    integer, intent(in) :: N_a,N_rho,N_t
    !Same number of gridpoints for agent 1 and agent 2 is assumed - this can be
    !generalized easily, but its difficult to see the benefit of doing it.
    integer :: j_a1,j_a2,j_rho,j_t !used for storing which points are to be used in interpolation.
    real(dp) :: t_a1,t_a2,t_rho,t_t ! weights in linear combinations t*x(j) + (1-t)*x(j+1)

    integer :: m_ind,i_ind !indices for cycling over elements of choice matrix.

    y = 0.0_dp !initilalize


    !This follows quite closely algorithm suggested in the Numerical recipes in Fortran book.
    !1) Find the points to be used in interpolation (for all variables, point with index j
    !is the point s.t. x lies between x(j) and x(j+1).
    call locate_gr(a1_gr,N_a,x(1),j_a1) !j_a1 will be the index of the grid point in a1_gr
    call locate_gr(a2_gr,N_a,x(2),j_a2)
    call locate_gr(rho_gr,N_rho,x(3),j_rho)

    !coefficients in linear combinations (value 0 is associated with the lower point)
    t_a1 = (a1_gr(j_a1+1)-x(1))/(a1_gr(j_a1+1)-a1_gr(j_a1))
    t_a2 = (a2_gr(j_a2+1)-x(2))/(a2_gr(j_a2+1)-a2_gr(j_a2))
    t_rho = (rho_gr(j_rho+1)-x(3))/(rho_gr(j_rho+1)-rho_gr(j_rho))

    !The interpolation for each of the elements of the consumption vector is done independently.

    !split The algorithm into 2 branches depending on whether we are also
    !interpolating over theta (static_environment = .false.) or not!
    if(N_t>1) then
        !4-dimensional interpolation
        call locate_gr(t_gr,N_t,x(4),j_t)
        t_t = (t_gr(j_t+1)-x(4))/(t_gr(j_t+1)-t_gr(j_t))

        y = C_pol_all(:,:,j_a1,j_a2,j_rho,j_t,s_ind)*t_a1*t_a2*t_rho*t_t+&
            C_pol_all(:,:,j_a1+1,j_a2,j_rho,j_t,s_ind)*(1-t_a1)*t_a2*t_rho*t_t+&
            C_pol_all(:,:,j_a1,j_a2+1,j_rho,j_t,s_ind)*t_a1*(1-t_a2)*t_rho*t_t+&
            C_pol_all(:,:,j_a1+1,j_a2+1,j_rho,j_t,s_ind)*(1-t_a1)*(1-t_a2)*t_rho*t_t+&
            C_pol_all(:,:,j_a1,j_a2,j_rho+1,j_t,s_ind)*t_a1*t_a2*(1-t_rho)*t_t+&
            C_pol_all(:,:,j_a1+1,j_a2,j_rho+1,j_t,s_ind)*(1-t_a1)*t_a2*(1-t_rho)*t_t+&
            C_pol_all(:,:,j_a1,j_a2+1,j_rho+1,j_t,s_ind)*t_a1*(1-t_a2)*(1-t_rho)*t_t+&
            C_pol_all(:,:,j_a1+1,j_a2+1,j_rho+1,j_t,s_ind)*(1-t_a1)*(1-t_a2)*(1-t_rho)*t_t+&
            C_pol_all(:,:,j_a1,j_a2,j_rho,j_t+1,s_ind)*t_a1*t_a2*t_rho*(1-t_t)+&
            C_pol_all(:,:,j_a1+1,j_a2,j_rho,j_t+1,s_ind)*(1-t_a1)*t_a2*t_rho*(1-t_t)+&
            C_pol_all(:,:,j_a1,j_a2+1,j_rho,j_t+1,s_ind)*t_a1*(1-t_a2)*t_rho*(1-t_t)+&
            C_pol_all(:,:,j_a1+1,j_a2+1,j_rho,j_t+1,s_ind)*(1-t_a1)*(1-t_a2)*t_rho*(1-t_t)+&
            C_pol_all(:,:,j_a1,j_a2,j_rho+1,j_t+1,s_ind)*t_a1*t_a2*(1-t_rho)*(1-t_t)+&
            C_pol_all(:,:,j_a1+1,j_a2,j_rho+1,j_t+1,s_ind)*(1-t_a1)*t_a2*(1-t_rho)*(1-t_t)+&
            C_pol_all(:,:,j_a1,j_a2+1,j_rho+1,j_t+1,s_ind)*t_a1*(1-t_a2)*(1-t_rho)*(1-t_t)+&
            C_pol_all(:,:,j_a1+1,j_a2+1,j_rho+1,j_t+1,s_ind)*(1-t_a1)*(1-t_a2)*(1-t_rho)*(1-t_t)

    else
        !3-dimensional interpolation (trilinear). In this case, the fourth index
        !is fixed, so we work with V(:,:,:,1,s_ind).

        y = C_pol_all(:,:,j_a1,j_a2,j_rho,1,s_ind)*t_a1*t_a2*t_rho + &
            C_pol_all(:,:,j_a1+1,j_a2,j_rho,1,s_ind)*(1-t_a1)*t_a2*t_rho +&
            C_pol_all(:,:,j_a1,j_a2+1,j_rho,1,s_ind)*t_a1*(1-t_a2)*t_rho +&
            C_pol_all(:,:,j_a1+1,j_a2+1,j_rho,1,s_ind)*(1-t_a1)*(1-t_a2)*t_rho +&
            C_pol_all(:,:,j_a1,j_a2,j_rho+1,1,s_ind)*t_a1*t_a2*(1-t_rho) +&
            C_pol_all(:,:,j_a1+1,j_a2,j_rho+1,1,s_ind)*(1-t_a1)*t_a2*(1-t_rho) +&
            C_pol_all(:,:,j_a1,j_a2+1,j_rho+1,1,s_ind)*t_a1*(1-t_a2)*(1-t_rho) +&
            C_pol_all(:,:,j_a1+1,j_a2+1,j_rho+1,1,s_ind)*(1-t_a1)*(1-t_a2)*(1-t_rho)


    end if


end subroutine interp_C

!Subroutine locate_gr is used to find the position of the point x at which we want to get V(x)
!in the grid. It is a slightly modified version of the code from Numerical Recipes in Fortran, 2nd edition.
!
!xx(1:n) is the grid, x is the point at which we want to get the value, j is such
!index that the point x lies between x(j) and x(j+1)
!
!The difference is that we now assume ascending order and we allow extrapolation.
!(The first change is for increase efficiency because this subroutine will be used
!very often, the second change is to help avoid errors due to numerical inaccuracies)
subroutine locate_gr(xx,n,x,j)
    integer, intent(in) :: n
    integer, intent(out) :: j
    real(dp), dimension(1:n), intent(in) :: xx
    real(dp), intent(in) :: x

    integer :: jl,jm,ju

    jl = 0
    ju = n+1

    do while(ju - jl > 1)
        jm = (ju+jl)/2 !mid point of bisection

        if(x>=xx(jm)) then
            jl = jm
        else
            ju = jm
        end if
    end do

    if(x<=xx(1)) then
        j=1
    else if(x>=xx(n)) then
        j=n-1
    else
        j=jl
    end if

end subroutine locate_gr


!subroutine LFFC (look for feasible choice) attempts to find a choice of consumption
!such that no constraint is violated, including the constraints imposed by grid range
!of a and rho. It does so for a given value of states (in the form used in generating the grid).

! Extension: If no such point found, return point where the violation was
!the lowest!!! For this I will need to properly finish implementation of the loss functions
!in check_constr subroutine.
subroutine LFFC(a_st,rho_st,t_st,g_st,g_st_ind,pars,b_min_mat,b_max_mat,l_max_mat,a_min_mat,a_max_mat,&
rho_min_vec,rho_max_vec,LFFC_choices,num_fc,WB_poss,c_guess_all_gp,fc_found,bcsf_loss,CTS_ind,C_pol_all,grids_tmp)
    Use nag_library, Only: e04jyf, nag_wp, e04jcf, e04jcp,e05jaf,e05jbf,E05JBK,E04ABA

    !Values of state variables as saved in the grid (hence the suffix st).
    real(dp), dimension(I_par), intent(in) :: a_st
    real(dp), dimension(I_par-1), intent(in) :: rho_st
    real(dp), intent(in) :: t_st,g_st
    integer, intent(in) :: g_st_ind
    type(par), intent(in) :: pars
    real(dp), dimension(M_par,I_par), intent(in) :: l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat
    real(dp), dimension(M_par,1), intent(in) :: rho_min_vec,rho_max_vec
    real(dp), dimension(:,:,:), intent(in) :: LFFC_choices
    integer, intent(in) :: num_fc
    integer, dimension(1+I_par*M_par,3) :: WB_poss !possibilities of what_binds
    type(c_guess_all), intent(out) :: c_guess_all_gp !output: all feasible choices at a grid point
    logical, intent(out) :: fc_found !true if a feasible choice found at this gridpoint
    type(grids_type), intent(in) :: grids_tmp

    real(dp), intent(out) :: bcsf_loss !minimum loss (that will either be zero if strictly feasible
    !point found or the loss at the minimum loss point that was found).

    integer, intent(in) :: CTS_ind !stage of CTS algorithm.
    real(dp), dimension(:,:,:,:,:,:,:), allocatable, intent(in) :: C_pol_all !policy function (from last stage of CTS algorithm)

    real(dp), dimension(1,I_par) :: theta_lvl
    real(dp), dimension(M_par,I_par) :: max_cons_mat
    real(dp), dimension(M_par,I_par) :: g_lvl_mat
    integer :: m_ind
    real(dp), dimension(M_par,M_par) :: P_onerow
    logical :: c_im_found

    integer :: WB_ind,wb_ind_max, fc_ind !what_binds index and index for feasible choice (used to loop through LFFC_choices)

    logical :: constr_ok
    real(dp) :: loss

    integer :: nf

    integer :: optim_subroutine = 1 !Temporarily place here,
    !later on put into parameter file. 1 means e04jcf (works on HPC),
    !2 means e04jyf

    !variables used for finding best choice (the least infeasible in the sense of having
    !the lowest loss) (bcsf stands for best choice so far)
    integer :: bcsf_wb_ind !what_binds index
    integer, dimension(3) :: bcsf_wb !what_binds for the best choice so far
    real(dp), dimension(M_par,I_par) :: bcsf_c_gp

    real(dp), dimension(M_par,I_par) :: l,a_prime,b_prime
    !real(dp), dimension(1,I_par) :: theta_prime

    !Local choice variables: gp stands for grid point.
    real(dp), dimension(M_par,I_par) :: c_gp

    integer :: fail_counter
    integer :: att_max !maximum number of attempts (stop when we find a strictly feasible choice)
    !if we do not find a feasible choice use the best of the att_max attempts in terms of loss (not adjusted - see
    !subroutine check_constr for loss adjustment meaning).
    !Need to experiment with this - large values might lead to too much time being spent in LFFC and it is not
    !clear how much a slightly better initial guess improves performance later in the program.
    integer :: att_counter !attempt counter (number of attempts tried).

    real(dp), dimension(ncvar_par) :: xc !choice variables (M_par x I_par matrix in row form)
    real(dp) :: fc !function value in maximization
    !The following two arrays are arrays of reals and integers which can be used to supply information
    !to the maximized function as an alternative to COMMON block.
    real(dp), dimension(1) :: ruser = 0.0_dp !Contains no information, only defined for purposes of using a NAG subroutine
    integer, dimension(3) :: iuser = 0 !will be used to pass what_binds to some subroutines.
    integer :: ifail = 0
    !Bounds:
    real(dp), dimension(M_par*I_par-1) :: bl = 0.0001_dp,bu
    real(dp) :: fmax
    !Variables used in 1-D optimization
    real(dp) :: E1,E2,A,B
    integer :: maxcal


    !declarations of copies of variables used in the common block
    real(dp), dimension(I_par) :: a_st_cp
    real(dp), dimension(I_par-1) :: rho_st_cp
    real(dp) :: t_st_cp,g_st_cp
    integer :: g_st_ind_cp
    real(dp), dimension(1,I_par) :: theta_lvl_cp
    real(dp), dimension(M_par,I_par) :: g_lvl_mat_cp
    real(dp), dimension(M_par,M_par) :: P_onerow_cp
    real(dp), dimension(M_par,I_par) :: l_max_mat_cp,b_min_mat_cp,b_max_mat_cp,a_min_mat_cp,a_max_mat_cp
    real(dp), dimension(M_par,1) :: rho_min_vec_cp,rho_max_vec_cp

    real(dp) :: loss_adj_min !minimized adjusted (to account for strictly satisfying some constraints) loss

    !The following common blocks are for passing data to loss_func (which needs to have a special interface
    !for use with NAG). We create local copies of the relevant variables. The reason for this is that
    !we pass some of these to LFFC so we might need to use the common block one level higher
    !instead of doing that, i.e., in subroutine get_c_guess_img. But this would not be very clean and would possibly
    !yield errors (no intent(in) checks). Also it will be easier to use this in TaxHetGrReg and other projects
    !if the adjustment happens in LFFC only and not on a higher level.
    !Creating a copy of the data is almost costless (small objects and the algorithm takes much more time in the iterative
    !minimization of loss_func.
    common /loss_func_bloc/ a_st_cp,theta_lvl_cp,rho_st_cp,t_st_cp,g_st_cp,g_st_ind_cp
    common /loss_func_bloc2/ g_lvl_mat_cp, P_onerow_cp,l_max_mat_cp,b_min_mat_cp,b_max_mat_cp&
    ,a_min_mat_cp,a_max_mat_cp,rho_min_vec_cp,rho_max_vec_cp

    optim_subroutine = pars%opt_subroutine

    !If the number of free variables in the optimization problem is 1, then
    !we can't use subroutine e04jcf. This will be the case only when we are solving
    !deterministic version of the model with two players.
    if(M_par*I_par-1 == 1) then
        optim_subroutine = 3
    end if

    !set number of attempts depending on which optimization routine we use
    att_max = pars%LFFC_attmax
    if(optim_subroutine == 3) att_max = 1 !the 1-d optimization subroutine does not depend on initial guess.


    !Start by computing maximum feasible consumption using max_feas_cons -> this will be the same
    !for all choices at this gridpoint and will be used to compute the actual choices.
    call prod_lvl(t_st,pars%theta_0,pars%xi_theta,pars%static_environment,theta_lvl)

    !Need government expenditure in levels for all shock realizations. Actually save it in a MxI matrix
    !where every row contains gov. exp. in level in state m
    do m_ind = 1,M_par
        g_lvl_mat(m_ind,:) = pars%S(1,m_ind)*trend_g(pars%mass,theta_lvl,pars%l_max,pars%k_gov)
    end do

    call max_feas_cons(pars%l_max,pars%mass,theta_lvl,g_lvl_mat,max_cons_mat)


    P_onerow = matmul(reshape([1.0_dp,1.0_dp],[M_par,1]),reshape(pars%P(g_st_ind,:),[1,M_par]))

    !Search for a feasible choice for all possible what_binds, that is, search for interior
    !solutions and all possible corner solutions.
    fc_found = .false. !false until at least one feasible choice found
    !initialize minimum loss
    bcsf_loss = 100000000000000000.0_dp


    !If CTS_ind > 1, we should have the policy function from the previous iteration (courser grid).
    !Use this to obtain initial guess by interpolation. If the initial guess is bad (in the sense of high loss),
    !then proceed further and perform the standard grid search for a feasible initial guess.
    if(CTS_ind>1 .and. pars%CTS_C_interp) then
        !Get the guess and save it into bcsf_c_gp (best choice so far in the sense of minimal loss)
        call interp_C(C_pol_all,[a_st(1),a_st(2),rho_st(1),t_st],g_st_ind,c_gp,grids_tmp%a1_gr&
            ,grids_tmp%a2_gr,grids_tmp%rho_gr,grids_tmp%t_gr,size(grids_tmp%a1_gr),size(grids_tmp%rho_gr),size(grids_tmp%t_gr))


        !The following computes loss. It's the same as in the part of the subroutine below so comments
        !were deleted.
        c_im_found = .false.
        WB_ind = 1 !interior solution
        call get_c_im(c_gp,WB_poss(WB_ind,:),a_st,rho_st,theta_lvl,g_lvl_mat(:,1)&
            ,g_st_ind,pars,P_onerow,c_im_found)


        !c_im_found will be .false. if the iterative algorithm couldn't find c from constraint
        !set loss very high
        if(c_im_found == .false.) then
            !c_gp(M_par,I_par) = 0.001_dp
            call fix_c(c_gp,rho_st,pars)
            go to 132
        end if

        call lab(c_gp,l,theta_lvl,pars%mass,g_lvl_mat(:,1),pars%Gamma_par)
        call ass(a_st,c_gp,l,pars%A_par,pars%B_par,pars%Gamma_par,pars%beta,P_onerow,a_prime)
        call get_b_prime(b_prime,a_prime,c_gp,pars)
        call check_constr(c_gp,l,b_prime,b_min_mat,b_max_mat,l_max_mat,theta_lvl,a_prime,a_min_mat&
            ,a_max_mat,WB_poss(WB_ind,:),rho_min_vec,rho_max_vec,constr_ok,loss,pars%A_par,fc_found,.false.)

        !Now check loss - if it's sufficiently small, consider the point feasible and
        132 if(loss < 0.05) then
                c_guess_all_gp%c_guess_all(1)%c_guess = c_gp
                c_guess_all_gp%c_guess_all(1)%what_binds = WB_poss(1,:)
                bcsf_loss = loss
                fc_found = .true. !Even slightly infeasible choices will be
                !considered feasible in this case.
                !If one feasible choice is enough (i.e., we try corner solutions only in cases
                !where we couldn't find a feasible choice for interior solution), return
                if(pars%one_fc_suff) return
        else !loss large, therefore the choice is not considered feasible. In that case
        !save if as best choice so far (including the loss) and proceed. If nothing better is found in the
        !grid search, then this choice may still be used as initial guess.

                bcsf_loss = loss
                bcsf_wb = WB_poss(1,:)
                bcsf_wb_ind = 1
                bcsf_c_gp = c_gp
        end if
    end if
    !End of the branch of program getting guess by interpolation from last CTS stage policy function


    !If we look for interior solutions only, only check the first WB possibility which corresponds to
    !interior solution
    if(pars%interior_solution_only) then
        wb_ind_max = 1
    else
        wb_ind_max = 1+I_par*M_par
    end if

    do WB_ind = 1,wb_ind_max

    !....Developing the improved LFFC algorithm.
    fail_counter = 0 !used for tracking how many times we were not able to find the missing element
    !of the cosumption vector.
    att_counter = 0 !initialize number of attempt.

    iuser = 0
    if(WB_ind > 1) then
        write(*,*) 'Need to generalize assignment of iuser in subroutine LFFC for WB_ind > 1 (corner solutions).'
    end if
    !(it should contain value of what_binds for the current WB_ind).

    !Before getting an initial guess of consumption and proceeding with the loss minimization, create
    !a copy of all the variables in the common blocks shared by subroutine loss_func.
    a_st_cp = a_st
    theta_lvl_cp = theta_lvl
    rho_st_cp = rho_st
    t_st_cp = t_st
    g_st_cp = g_st
    g_st_ind_cp = g_st_ind
    g_lvl_mat_cp = g_lvl_mat
    P_onerow_cp = P_onerow
    l_max_mat_cp = l_max_mat
    b_min_mat_cp = b_min_mat
    b_max_mat_cp = b_max_mat
    a_min_mat_cp = a_min_mat
    a_max_mat_cp = a_max_mat
    rho_min_vec_cp = rho_min_vec
    rho_max_vec_cp = rho_max_vec

    !get upper bounds bu for minimization of loss over consumption.
    call max_feas_cons(pars%l_max,pars%mass,theta_lvl,g_lvl_mat,max_cons_mat)
    !Use this to generate the vector of upper bounds...
    call c_inc_mtv(max_cons_mat,bu)



    !Generate random initial choice of consumption.
    150 call random_number(c_gp)

    att_counter = att_counter + 1 !increase the counter of attempts made

    c_gp = c_gp*max_cons_mat*0.3_dp !the 0.3 scaling is used because for values close to max_cons_mat, those
    !choice will not be feasible anyway so it's a waste of time to try these
    !c_gp contains the actual consumption choices by agents (corresponding to c's in the paper).


    call get_c_im(c_gp,WB_poss(WB_ind,:),a_st,rho_st,theta_lvl,g_lvl_mat(:,1)&
    ,g_st_ind,pars,P_onerow,c_im_found)
    if(c_im_found == .false. .or. .true.) call fix_c(c_gp,rho_st,pars)

    !________________________End of fixing initial guess___________________

    !Transform the initial guess into a vector (from a matrix) so that we can use it in loss_func and it's variants.
    call c_inc_mtv(c_gp,xc)

    ifail = 1 !silent exit (we do not want the program to stop when soft errors happen - for example sometimes
    !it is bound to happen that e04jcf does not find the optimum within the limit of iterations).

    select case(optim_subroutine)
        case(1) !Quadratic approximation method

            call e04jcf(loss_func2,M_par*I_par-1,ncvar_par*2+1,xc,bl,bu,0.1_dp,&
                pars%rhoend,e04jcp,pars%LFFC_maxiter,loss_adj_min,nf,iuser,ruser,ifail)

        case(2) !Newton-like method
            !call e04jyf(M_par*I_par-1,0,CPC_ret,bl,bu,xc,fmax,iw,liw,w,lw,iuser,ruser,ifail)
            if(this_image()==1) then
                write(*,*) 'Error: need to generalize LFFC for optim_subroutine = 2'
                error stop
            end if

        case(3) !1-D optimization (only used when M_par = 1, I_par = 2)
            A = bl(1)
            B = bu(1)
            e1 = pars%e1
            e2 = pars%e2
            maxcal = pars%LFFC_maxiter

            call E04ABA(loss_func_1d,E1,E2,A,B,MAXCAL,XC(1),loss_adj_min,IUSER,RUSER,IFAIL)

        case default
            write(*,*) 'Error in subroutine VFI_solve_gp - wrong value of optim_subroutine'
            error stop
    end select

    !Reshape the result back into incomplete consumption matrix form
    call c_inc_vtm(xc,c_gp)

    !Now compute the loss of the point (unadjusted) to evaluate how good the point is in terms of feasibility

    !We need to get the other variables. First of all get the last element of the
    !consumption matrix. This depends on what_binds.
    c_im_found = .false.
    call get_c_im(c_gp,WB_poss(WB_ind,:),a_st,rho_st,theta_lvl,g_lvl_mat(:,1)&
    ,g_st_ind,pars,P_onerow,c_im_found)

    !In deterministic case and interior solution, we are guaranteed to find a consumption plan that satistifies the non-
    !negativity constraints for all agents.
    if(c_im_found == .false.) then

        !Try a few times with different initial guess. If it still
        !does not work, just use zero for c_gp(M_par,I_par)
        fail_counter = fail_counter + 1
        if(fail_counter < 10) then
            go to 150 !Try again with a different randomly generated point.
        else
            !Everything else failed so use subroutine fix_c to get an initial guess
            !which may be far from optimum and may violate other constraints, but at least
            !all elements of consumption plan are non-negative.
            call fix_c(c_gp,rho_st,pars)
        end if
    end if

    !Now after we got initial choice just get next-period states and try computing the loss (generalized,
    !i.e., negative for strictly feasible choices). After that works the computing of loss will be
    !done by a function of c_gp (just like CPC_ret interface but returns loss, not current plus continuation
    !return) - so the things like getting next-period states and computing the loss will be moved into this
    !function.
    !Get labour supply implied by the consumption
    call lab(c_gp,l,theta_lvl,pars%mass,g_lvl_mat(:,1),pars%Gamma_par)
    !next-period states
    call ass(a_st,c_gp,l,pars%A_par,pars%B_par,pars%Gamma_par,pars%beta,P_onerow,a_prime)
    call get_b_prime(b_prime,a_prime,c_gp,pars)

    !Compute the unadjusted loss (so the last argument of check_constr is false). This is used for ranking
    !the points because feasibility has priority. In any case, whether we should the adjusted or unadjusted
    !loss shouldn't make much difference (they are closely related and in the case where ranking reversed,
    !we would still expect the points to be close and the maximization in later stages of the algorithm should
    !thus converge to a similar point even in the presence of nonconvex feasible sets.).
    call check_constr(c_gp,l,b_prime,b_min_mat,b_max_mat,l_max_mat,theta_lvl,a_prime,a_min_mat&
    ,a_max_mat,WB_poss(WB_ind,:),rho_min_vec,rho_max_vec,constr_ok,loss,pars%A_par,fc_found,.false.)

    !If the choice is feasible, save it and go to 200 (try finding feasible choices for next value
    !of what_binds (if indeed we care about other solutions than interior ones)
    if(constr_ok) then
        fc_found = .true. !we found at least one feasible choice - this means that
        !from now on we only care aboud strictly feasible choices (for other values of what_binds)
        !Save the result into c_guess_all_gp which contains all feasible guesses
        !at a gridpoint (or one as little infeasible guess as possible).
        c_guess_all_gp%c_guess_all(wb_ind)%c_guess = c_gp
        c_guess_all_gp%c_guess_all(wb_ind)%what_binds = WB_poss(WB_ind,:)
        bcsf_loss = loss !this should be zero in case of strictly feasible points.

        !If we were only looking for interior solution or if we are satisfied with
        !one feasible choice per gridpoint, return
        if(pars%one_fc_suff) return
        !go to next what_binds possibility
        go to 200

    else !constr_ok == .false., so choice not strictly feasible.

        !bcsf_loss will be 0 if a feasible choice was found so this thing is performed
        !only if we did not find a strictly feasible solution! We're already in the branch
        !of program where we know that the choice is not strictly feasible, so strictly
        !feasible choices are never ignored even if we found one before for different what_binds
        !unless we have pars%one_fc_suff == .true.
        if(loss < bcsf_loss) then
            bcsf_loss = loss
            bcsf_wb = WB_poss(WB_ind,:)
            bcsf_wb_ind = wb_ind
            bcsf_c_gp = c_gp
        end if

    !If the number of attempts is less than maximum number of attempts, try again with a different random initial
    !guess
    if(att_counter < att_max) go to 150

    end if

200 end do



end subroutine LFFC

!subroutine LFFC_msg displays results of looking for feasible choice (across all images)
subroutine LFFC_msg(share_feas_img,avg_loss_img,min_ind_img,max_ind_img,spec_img,V_num_el,folder_name,initial_time)
    real(dp), codimension[*] :: share_feas_img, avg_loss_img
    integer, intent(in) :: min_ind_img,max_ind_img,spec_img,V_num_el
    character*(*) :: folder_name
    real(dp), intent(in) :: initial_time

    real(dp) :: share_feas,avg_loss
    integer :: img_index

    character(256) :: tmp_string,tmp_string2,elapsed_time_string

!weigh the feasible shares by the number of gridpoints that the image serves
share_feas_img = share_feas_img*(max_ind_img-min_ind_img + 1.0_dp)
avg_loss_img = avg_loss_img*(max_ind_img-min_ind_img + 1.0_dp)

!(EFE) could avoid this sync all? Shouldn't really matter though because this thing immediately follows
!another sync all so pretty much no delay.

sync all !make sure that all images weighed their share
if(this_image() == spec_img) then
    share_feas = 0.0_dp
    avg_loss = 0.0_dp
    do img_index = 1,num_images()
        share_feas = share_feas + share_feas_img[img_index]
        avg_loss = avg_loss + avg_loss_img[img_index]
    end do
    share_feas = share_feas/real(V_num_el)
    avg_loss = avg_loss/real(V_num_el)

    call double_output(folder_name,'Finished looking for feasible choices.')
    write(tmp_string2,'(F6.4)') share_feas
    tmp_string = 'Share of gridpoints with feasible choice = '//trim(tmp_string2)
    call double_output(folder_name,trim(tmp_string))
    write(tmp_string2,'(F15.3)') avg_loss
    tmp_string = 'Average loss function value = '//trim(tmp_string2)
    call double_output(folder_name,trim(tmp_string))


    call runtime_report(initial_time,.true.,elapsed_time_string)
    call double_output(folder_name,trim(elapsed_time_string))
    call double_output(folder_name,' ')
end if

end subroutine LFFC_msg


!Subroutine LFFC_prep_choices prepares choices of consumption which will be tried
!in subroutine LFFC (which looks for feasible choices at a given grid point).
!
!The form in which the choices are stored is share of maximum feasible consumption in
!a state by an agent. The reason to do this at the outset and than use the same choices
!all over again is reducing runtime. This subroutine could be parallelized but it is not worth it
!becasue it is not a bottleneck - for such high value of LFFC_max_ind that this subroutine
!would run long, the looking for feasible choices itself would take prohibitively long time.
!
!The choices returned in LFFC_choices are already normalized by mass of agents, so to get the actual
!consumption choice we only need to multiply the maximum feasible consumption of all agents by these
!numbers.
subroutine LFFC_prep_choices(pars,LFFC_choices,num_fc)
    type(par), intent(in) :: pars
    real(dp), dimension(:,:,:) :: LFFC_choices !First 2 indices will be for MxI choice matrix,
    !last index corresponds to 'big' index
    integer, intent(out) :: num_fc !number of feasible choices found

    real(dp), dimension(pars%LFFC_max_ind) :: share_sp !share space

    integer :: num_cv,num_tmp !number of choice variables, and a temporary constant which is useful to pre-compute
    integer :: k,m,R
    real(dp) :: share_sum !sum of shares on maximum feasible consumption to group (must be <= 1)
    integer, dimension(I_par*M_par-1) :: indices !contains indices for all the choice variables
    !which we will retrieve from the one big index k

    real(dp), dimension(I_par*M_par-1,1) :: LFFC_choice_one !One particular combination of choice variables

    !First find a space of shares which will be then used to look for the feasible choices.
    !If this space is changed, we would then change the order in which the
    !algorithm is looking for feasible choices, which may be important given the
    !non-convexity of the constraint set. I leave it as a generalization - we can also have
    !different space of shares for different choice variables.

    !The space needs not be linear - in fact it might be more likely to find feasible choices
    !for relatively low shares of maximum consumption.

    !all parameters governing the space of shares used for choices generation are in parameter file
    call gen_equi_grid(share_sp,pars%LFFC_min_share,pars%LFFC_max_share)

    !pre-compute some values used lots of times
    num_cv = I_par*M_par-1
    num_tmp = I_par*(M_par-2)+1 !The last index at which summation should start in the 'big vector'

    num_fc = 0 !initialize number of feasible choice found
    do k = 1,pars%LFFC_max_ind**num_cv
        R=k !initialize remainder
        do m=1,num_cv
            indices(m) = (R-1)/(pars%LFFC_max_ind**(num_cv-m)) + 1
            R=R-(indices(m)-1)*(pars%LFFC_max_ind**(num_cv-m))
        end do
        !Now we have the indices so we can compute the choice, but save it only
        !if feasibility constraint is not violated. If that is the case,
        !increment num_fc by 1 and use that so save results, so that num_fc
        !at the end will be the number of feasible choices found.
        LFFC_choice_one(:,1) = share_sp(indices)

        do m=1,num_tmp,I_par
            share_sum = sum(LFFC_choice_one(m:m+I_par-1,1))
            !Here we are not concerned with he mass of the agents, just the share going to each group
            !If at any point we find summation of share bigger than one, do not save this choice.

            if(share_sum > 1.0_dp) goto 100
        end do

        !Saving of results goes here - if sum of shares bigger than 1 then this part is skipped using goto
        num_fc = num_fc + 1

        !Need to convert the choice into a matrix. Just be careful about order of elements.

        !The missing element will become a zero.
        LFFC_choices(:,:,num_fc) = transpose(reshape(LFFC_choice_one(:,1),[I_par,M_par],[0.0_dp]))

100 end do
!debug:

!write(*,*) real(num_fc)/real(pars%LFFC_max_ind**num_cv)

!What we worked with until now are shares of consumption of all (potentially) available resources
!consumed by a group of agents (so that c_i*pi_i = something). It is more convenient to save these
!choices already normalized by mass of agents, so we can get the actual consumption choices
!simply by multiplying the saved share by the maximum consumption, and we won't have to divide by mass
!every time.
!Do this for every player
do m=1,I_par
    LFFC_choices(:,m,:) = LFFC_choices(:,m,:)/pars%mass(1,m)
end do

end subroutine LFFC_prep_choices

!subroutine load_N loads the number of gridpoints. This is then used to allocate arrays for
!loading the other variables. M_par_loc is loaded only to check that the value
!function in the file was generated for the same value of M_par, otherwise it makes no sense to use
!the value function loaded from file.
subroutine load_N(N_a,N_rho,N_t,input_folder,folder_name)
    integer, intent(out) :: N_a,N_rho,N_t

    character*(*), intent(in) :: input_folder,folder_name !folder_name is the output folder (in which log is saved)
    character*256 :: shell_cmd,file_path

    integer :: M_par_loc

    logical :: file_exists


    !First check whether the results are in the folder. Of course the program would crash if we tried to read from
    !non-existent files but this will help locating the problem.
    inquire(file=('results/'//trim(input_folder)//'/N_a.out'), exist=file_exists)
    if (.not. file_exists) then
        call double_output(folder_name,'Error: File N_a.out in input_folder not found.')
        call double_output(folder_name,'Change input_folder or load_initial_guess.')
        error stop
    end if

    inquire(file=('results/'//trim(input_folder)//'/N_rho.out'), exist=file_exists)
    if (.not. file_exists) then
        call double_output(folder_name,'Error: File N_rho.out in input_folder not found.')
        call double_output(folder_name,'Change input_folder or load_initial_guess.')
        error stop
    end if

    inquire(file=('results/'//trim(input_folder)//'/N_t.out'), exist=file_exists)
    if (.not. file_exists) then
        call double_output(folder_name,'Error: File N_t.out in input_folder not found.')
        call double_output(folder_name,'Change input_folder or load_initial_guess.')
        error stop
    end if

    inquire(file=('results/'//trim(input_folder)//'/M_par_loc.out'), exist=file_exists)
    if (.not. file_exists) then
        call double_output(folder_name,'Error: File M_par_loc.out in input_folder not found.')
        call double_output(folder_name,'Change input_folder or load_initial_guess.')
        error stop
    end if

    !load the variables from files
    file_path = 'results/'//trim(input_folder)//'/N_a.out'
    open(unit=20, file = file_path, status = 'old')
    read(unit=20, fmt=*)  N_a
    close(unit = 20)

    file_path = 'results/'//trim(input_folder)//'/N_rho.out'
    open(unit=20, file = file_path, status = 'old')
    read(unit=20, fmt=*)  N_rho
    close(unit = 20)

    file_path = 'results/'//trim(input_folder)//'/N_t.out'
    open(unit=20, file = file_path, status = 'old')
    read(unit=20, fmt=*)  N_t
    close(unit = 20)

    file_path = 'results/'//trim(input_folder)//'/M_par_loc.out'
    open(unit=20, file = file_path, status = 'old')
    read(unit=20, fmt=*)  M_par_loc
    close(unit = 20)

    if(M_par_loc /= M_par) then
        write(*,*) 'Error: M_par_loc does not equal M_par.'
        write(*,*) 'Value function computed for different M_par is useless as initial guess.'
        error stop
    end if

end subroutine load_N

!subroutine load_V loads the value function and grids. Interpolation is used. Due to numerical errors,
!it is easiest to use interpolation every time even if the grid were identical.
subroutine load_V(V_vect,a1_gr,a2_gr,rho_gr,t_gr,input_folder,folder_name)
    real(dp), dimension(:), intent(out) :: V_vect
    real(dp), dimension(:), intent(out) :: a1_gr,a2_gr,rho_gr,t_gr

    character*(*), intent(in) :: input_folder,folder_name !folder_name is the output folder (in which log is saved)
    character*256 :: file_path

    logical :: file_exists

    !Check that all the grid files and the value function file exist.
    inquire(file=('results/'//trim(input_folder)//'/A1_gr.out'), exist=file_exists)
    if (.not. file_exists) then
        call double_output(folder_name,'Error: File A1_gr.out in input_folder not found.')
        call double_output(folder_name,'Change input_folder or load_initial_guess.')
        error stop
    end if
    inquire(file=('results/'//trim(input_folder)//'/A2_gr.out'), exist=file_exists)
    if (.not. file_exists) then
        call double_output(folder_name,'Error: File A2_gr.out in input_folder not found.')
        call double_output(folder_name,'Change input_folder or load_initial_guess.')
        error stop
    end if
    inquire(file=('results/'//trim(input_folder)//'/rho_gr.out'), exist=file_exists)
    if (.not. file_exists) then
        call double_output(folder_name,'Error: File rho_gr.out in input_folder not found.')
        call double_output(folder_name,'Change input_folder or load_initial_guess.')
        error stop
    end if
    inquire(file=('results/'//trim(input_folder)//'/t_gr.out'), exist=file_exists)
    if (.not. file_exists) then
        call double_output(folder_name,'Error: File t_gr.out in input_folder not found.')
        call double_output(folder_name,'Change input_folder or load_initial_guess.')
        error stop
    end if
    inquire(file=('results/'//trim(input_folder)//'/V_vect.out'), exist=file_exists)
    if (.not. file_exists) then
        call double_output(folder_name,'Error: File V_vect.out in input_folder not found.')
        call double_output(folder_name,'Change input_folder or load_initial_guess.')
        error stop
    end if

    !Reading all the variables
    file_path = 'results/'//trim(input_folder)//'/V_vect.out'
    open(unit=20, file = file_path, status = 'old', form='binary')
    read(unit=20)  V_vect
    close(unit = 20)

    file_path = 'results/'//trim(input_folder)//'/A1_gr.out'
    open(unit=20, file = file_path, status = 'old')
    read(unit=20, fmt=*)  A1_gr
    close(unit = 20)

    file_path = 'results/'//trim(input_folder)//'/A2_gr.out'
    open(unit=20, file = file_path, status = 'old')
    read(unit=20, fmt=*)  A2_gr
    close(unit = 20)

    file_path = 'results/'//trim(input_folder)//'/rho_gr.out'
    open(unit=20, file = file_path, status = 'old')
    read(unit=20, fmt=*)  rho_gr
    close(unit = 20)

    file_path = 'results/'//trim(input_folder)//'/t_gr.out'
    open(unit=20, file = file_path, status = 'old')
    read(unit=20, fmt=*)  t_gr
    close(unit = 20)

end subroutine load_V

!subroutine load_guess loads the initial guess of consumption allocation (which is
!assumed to have been saved by subroutine C_guess).
subroutine load_guess(c_guess_c,c_guess_wb,input_folder,folder_name)
    real(dp), dimension(:,:,:), intent(out) :: c_guess_c
    integer, dimension(:,:), intent(out) :: c_guess_wb

    character*(*), intent(in) :: input_folder,folder_name !folder_name is the output folder (in which log is saved)
    character*256 :: file_path

    !no checking whether the files exist (it was done before calling this subroutine).

    !Reading all the variables
    file_path = 'results/'//trim(input_folder)//'/c_guess_c.out'
    open(unit=20, file = file_path, status = 'old', form='binary')
    read(unit=20)  c_guess_c
    close(unit = 20)

    file_path = 'results/'//trim(input_folder)//'/c_guess_wb.out'
    open(unit=20, file = file_path, status = 'old', form='binary')
    read(unit=20)  c_guess_wb
    close(unit = 20)

end subroutine load_guess


!The following subroutine computes maximum feasible consumption going to a group
!of agents for given productivity and government expenditure in levels. To get
!the maximum feasible consumption by an agent, this still needs to be divided by
!the mass of the agent!
subroutine max_feas_cons(l_max,mass,theta_lvl,g_lvl,max_cons)
    real(dp), dimension(1,I_par), intent(in) :: l_max,theta_lvl,mass
    real(dp), dimension(M_par,I_par) ,intent(in) :: g_lvl

    real(dp),dimension(M_par,I_par), intent(out) :: max_cons

    integer :: m_ind

    !Reminder: This is maximum consumption by any agent, it's the same for both agents
    !and depends on state only.
    do m_ind = 1,M_par
        max_cons(m_ind,:) = sum(mass*theta_lvl*l_max) - g_lvl(m_ind,1)
    end do

end subroutine max_feas_cons


!Subroutine prep_bounds converts some of the bounds which are not dependent
!on state into state-by-state matrices of bounds. These are then used for checking
!constraints more efficiently than using loops or if clauses.
subroutine prep_bounds(l_max,b_min,b_max,a_min,a_max,rho_min,rho_max,l_max_mat,b_min_mat,b_max_mat,&
a_min_mat,a_max_mat,rho_min_vec,rho_max_vec)
    real(dp), dimension(1,I_par), intent(in) :: l_max,b_min,b_max,a_min,a_max
    real(dp), intent(in) :: rho_min,rho_max

    real(dp), dimension(M_par,I_par), intent(out) :: l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat
    real(dp), dimension(M_par,1), intent(out) :: rho_min_vec,rho_max_vec

    integer :: i

    do i=1,I_par
        l_max_mat(:,i) = l_max(1,i)
    end do
    do i=1,I_par
        b_min_mat(:,i) = b_min(1,i)
    end do
    do i=1,I_par
        b_max_mat(:,i) = b_max(1,i)
    end do
    do i=1,I_par
        a_min_mat(:,i) = a_min(1,i)
    end do
    do i=1,I_par
        a_max_mat(:,i) = a_max(1,i)
    end do

    rho_min_vec(:,1) = rho_min
    rho_max_vec(:,1) = rho_max

end subroutine prep_bounds

!subroutine prod_lvl computes productivities of all agents in levels from the
!time state variable. In time this will be generalized to allow arbitrary functions
!or lists of values directly - right now this assumes constant growth as the
!first paper.
subroutine prod_lvl(t_st,theta_0,xi_theta,static_environment,theta_lvl)
    real(dp), intent(in) :: t_st
    real(dp), dimension(1,I_par), intent(in) :: theta_0,xi_theta
    logical, intent(in) :: static_environment

    real(dp), dimension(1,I_par), intent(out) :: theta_lvl

    real(dp) :: t !time period

    if(static_environment == .true.) then
        theta_lvl = theta_0
        return
    end if

    theta_lvl = theta_0*((xi_theta+1.0_dp)**t_st)
end subroutine prod_lvl


function rel_prod_prime(rel_growth,rel_prod,rel_prod_min,static_environment)
    real(dp) :: rel_prod_prime
    real(dp), intent(in) :: rel_growth,rel_prod,rel_prod_min !relative growth, relative productivity,
    !and the lowest permitted relative productivity (this will be on the grid)
    logical, intent(in) :: static_environment

    !The case of static_environment == .true. is covered by the following formula because
    !this function never returns productivity lower than rel_prod_min
    !which will be the current state

    !However, due to possible numerical error and resulting bugs in the case of static_productivity,
    !return the exact value if that is that case.

    if(static_environment == .true.) then
        rel_prod_prime = rel_prod
        return
    end if

    rel_prod_prime = maxval([rel_growth * rel_prod,rel_prod_min])

end function


!Subroutine save_V saves the value function and the grid for which it was generated.
!It also saves the number of grid points. This is useful when we load the
!value function from file. M_par is saved just for checking purposes - in practice
!we will never want to use a value function computed for different M.
subroutine save_V(V_vect,grids,N_a,N_rho,N_t,M_par_loc,folder_name)
    real(dp), dimension(:), intent(in) :: V_vect
    type(grids_type), intent(in) :: grids
    integer, intent(in) :: N_a,N_rho,N_t,M_par_loc !M_par_loc is just a local copy

    character*(*), intent(in) :: folder_name
    character*256 :: shell_cmd,file_path

    !Saving all the variables
    file_path = 'results/'//trim(folder_name)//'/V_vect.out'
    open(unit=20, file = file_path, status = 'replace', form='binary') !experiment with form = 'binary'
    write(unit=20)  V_vect
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/A1_gr.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  grids%A1_gr
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/A2_gr.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  grids%A2_gr
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/rho_gr.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  grids%rho_gr
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/t_gr.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  grids%t_gr
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/N_a.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  N_a
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/N_rho.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  N_rho
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/N_t.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  N_t
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/M_par_loc.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  M_par_loc
    close(unit = 20)


end subroutine save_V


!Subroutine save_guess saves one guess for every gridpoint, i.e., the
!consumption allocation and the what_binds. These arrays are saved
!in two separate arrays
subroutine save_guess(c_guess_c,c_guess_wb,folder_name)
    real(dp), dimension(:,:,:), intent(in) :: c_guess_c
    integer, dimension(:,:), intent(in) :: c_guess_wb


    character*(*), intent(in) :: folder_name
    character*256 :: shell_cmd,file_path,tmp_string

    file_path = 'results/'//trim(folder_name)//'/c_guess_c.out'
    open(unit=20, file = file_path, status = 'replace', form='binary')
    write(unit=20)  c_guess_c
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/c_guess_wb.out'
    open(unit=20, file = file_path, status = 'replace', form='binary')
    write(unit=20)  c_guess_wb
    close(unit = 20)

end subroutine save_guess



!SWF(c,l,s_current) is expected social welfare, given a matrix of state-contingent
!consumption and labour.
!Each row of c and l corresponds to a state, each column to an agent. s_ind_st is the
!index of last-period shock realization (state). mass and alpha are a vector of
!mass of agents and Pareto weights in the social welfare function.
function SWF(c,l,s_ind_st,P,mass,alpha,A_par,B_par,Gamma_par)
    real(dp) :: SWF

    real(dp), dimension(M_par,I_par), intent(in) :: c,l
    real(dp), dimension(M_par,M_par), intent(in) :: P
    integer, intent(in) :: s_ind_st
    real(dp), dimension(1,2), intent(in) :: mass,alpha
    real(dp), intent(in) :: A_par,B_par,Gamma_par !parameters of the utility function

    real(dp), dimension(M_par,I_par) :: u !Utility of all agents in all states.

    call util(c,l,U,A_par,B_par,Gamma_par) !need to use mod_util for this

    SWF = maxval(matmul(matmul(P(s_ind_st:s_ind_st,1:M_par),U),transpose(mass*alpha)))
    !The maxval is there to convert array of dimension (1,1) to scalar. Perhaps this is not the
    !most efficient way to do this but I don't know how else to do it in a single statement,
    !i.e., without saving the result as (1,1) array and then accesing this array.
    !
    !This function is called potentially hundreds of millions of time per execution so maximum efficiency is important here.
end function SWF





!Function trend_g computes trend government expenditure.
function trend_g(mass,theta_lvl,l_max,k_gov)
    real(dp) :: trend_g
    real(dp), dimension(1,I_par), intent(in) :: mass,theta_lvl, l_max
    real(dp), intent(in) :: k_gov

    trend_g = k_gov*sum(mass*theta_lvl*l_max)

end function trend_g


!Subroutine V_from_V gets the value of Value function at all gridpoints
!from a different value function using interpolation. The input value function is denoted
!V_old (with associted grids grids_old), the output is V_new with grids grids_new (the grids
!are also an input).
!
!The value functions are in multi-dimensional array format, not in row format,
!so this function will need to be generalized by switching to row format when
!we want to solve the model for I>2. This would make it less efficient though
!(or we would also have to generalize the interpolation function) and is
!not worth the effort at this time.
subroutine V_from_V(V_old,grids_old,grids_new,V_new,folder_name)
    real(dp), dimension(:,:,:,:,:), intent(in) :: V_old !5 states if I=2
    type(grids_type), intent(in) :: grids_old,grids_new

    real(dp), dimension(:,:,:,:,:), intent(out) :: V_new

    character*(*), intent(in) :: folder_name

    integer, dimension(5) :: shape_old,shape_new,shape_grid_old,shape_grid_new

    integer :: a1_ind,a2_ind,rho_ind,t_ind,g_ind

!Perform checks that the grids and the value functions are compatible in the sense of dimensionality.

!(EFE) can make check that if the grids are identical, V_old will simply replace V_new
!This will be very race case though, only happening if we don't use CTS algorithm.
!can save a bit of time with really large grids

    shape_old = shape(V_old)
    shape_new = shape(V_new)

    shape_grid_old(1) = size(grids_old%a1_gr)
    shape_grid_old(2) = size(grids_old%a2_gr)
    shape_grid_old(3) = size(grids_old%rho_gr)
    shape_grid_old(4) = size(grids_old%t_gr)
    !(grids for shock realization not checked). They should always
    !be the same and we only check the whether the value functions
    !have the same M_par. Also, no interpolation happens over shock
    !realizations - we just assume that shock space is the same here.

    shape_grid_new(1) = size(grids_new%a1_gr)
    shape_grid_new(2) = size(grids_new%a2_gr)
    shape_grid_new(3) = size(grids_new%rho_gr)
    shape_grid_new(4) = size(grids_new%t_gr)

    if(any(shape_old(1:4)/=shape_grid_old(1:4))) then
        call double_output(folder_name, 'Error: grids_old in subroutine V_from_V incompatible with V_old.')
        error stop
    end if
    if(any(shape_new(1:4)/=shape_grid_new(1:4))) then
        call double_output(folder_name, 'Error: grids_new in subroutine V_from_V incompatible with V_new.')
        error stop
    end if

    !If any of the gridpoints has different M_par, give an error. These should never happen but check
    !just in case, because a bug here could be difficult to detect.
    if(shape_old(5) /= M_par) then
        call double_output(folder_name, 'Error: V_old in V_from_V computed for different M_par than the current one.')
        error stop
    end if
    !If any of the gridpoints has different M_par, give an error.
    if(shape_new(5) /= M_par) then
        call double_output(folder_name, 'Error: V_new in V_from_V computed for different M_par than the current one.')
        error stop
    end if


    !ALSO check whether the ranges of the new grid exceed those of the old grid. If yes,
    !extrapolation is used. For small differences this is fine but a warning should be displayed.
    !(do not assume monotonicy in grids)
    if((minval(grids_new%a1_gr)<minval(grids_old%a1_gr)) .or. &
    (minval(grids_new%a2_gr)<minval(grids_old%a2_gr)) .or. &
    (minval(grids_new%rho_gr)<minval(grids_old%rho_gr)) .or. &
    (minval(grids_new%t_gr)<minval(grids_old%t_gr)) .or. &
    (maxval(grids_new%a1_gr)>maxval(grids_old%a1_gr)) .or. &
    (maxval(grids_new%a2_gr)>maxval(grids_old%a2_gr)) .or. &
    (maxval(grids_new%rho_gr)>maxval(grids_old%rho_gr)) .or. &
    (maxval(grids_new%t_gr)>maxval(grids_old%t_gr))) then
        call double_output(folder_name,'Warning: in V_from_V, bounds of grid_new exceed bounds of grid_old.')
        call double_output(folder_name,'Extrapolation will be used, dangerous far from bounds.')
    end if

    !Later on, can use some NAG functions. So far use the interpolation functions which I wrote.

    !When we generalize to I>2, there will no longer be exactly 5 state variables. In this case,
    !we will have only one main index here and we will use ind_unfold_all to cycle over the new value function.
    !(Or maybe just generate is as when ind_unfold_all was generated on the go because we saved
    !only indices for the current image to save memory).
    !Or generate the unfolded indices for all images here at the start and then cycle.


    do a1_ind = 1,shape_grid_new(1)
        do a2_ind = 1,shape_grid_new(2)
            do rho_ind = 1,shape_grid_new(3)
                do t_ind = 1,shape_grid_new(4)
                    do g_ind = 1,M_par
                    !get the interpolated value and save it in the appropriate element of V_new
                    call interp_V(V_old,[grids_new%a1_gr(a1_ind),grids_new%a2_gr(a2_ind),grids_new%rho_gr(rho_ind),&
                        grids_new%t_gr(t_ind)],g_ind,V_new(a1_ind,a2_ind,rho_ind,t_ind,g_ind)&
                        ,grids_old%a1_gr,grids_old%a2_gr,grids_old%rho_gr,grids_old%t_gr,&
                        shape_grid_old(1),shape_grid_old(3),shape_grid_old(4))
                end do
            end do
        end do
    end do
end do
end subroutine V_from_V

!Subroutine VFI_solve_gp solves the value function iteration maximization step at a single
!grid point. This subroutine needs NAG libraries.
subroutine VFI_solve_gp(ind_unfold,grids,pars,N_a,N_rho,N_theta,c_guess_all_gp,V_new_gp,skip_max,&
C_pol_gp,glob_opt,loss_gp,optim_wb_gp)
    Use nag_library, Only: e04jyf, nag_wp, e04jcf, e04jcp,e05jaf,e05jbf,E05JBK,E04ABA
    !USE nag_library

    integer, dimension(5), intent(in) :: ind_unfold !vector of indices (this assumes I=2, generalization will be needed)
    type(grids_type), intent(in) :: grids
    type(par), intent(in) :: pars
    integer, intent(in) :: N_a,N_rho,N_theta !number of gridpoints (no need to compute them here again)
    type(c_guess_all), intent(inout) :: c_guess_all_gp !contains the initial guess of choice variables at this gp
    !It will be used in maximization and the new optimal choice will be returned in its place
    real(dp), intent(out) :: V_new_gp !New optimum value at the grid point.
    logical, intent(in) :: skip_max !If true, then the maximization step will be skipped and the
    !function returns the same initial guess as was passed to it, and the value obtained
    !using this initial guess.
    logical, intent(in) :: glob_opt !If this is true, a global maximization subroutine will be used.
    !This is probably a good idea to do in the first iteration when the initial gues may be quite wrong
    !(And due to non-convexity of the feasible set we might end up at a global maximum far
    !from the actual global maximum). However, this is very slow so it should really
    !be done only in the first iteration (and then maybe once in a while as a robustness check).

    real(dp), dimension(M_par,I_par), intent(out) :: C_pol_gp !optimal choice at a grid point

    !integer, intent(out) :: failcount !serves for debugging purposes,
    !will be 0 on exit if a maximization routine finished without any issues,
    !if will be 1 if something bad happened (in the sense of a NAG subroutine
    !returning a serious error).

    !Pointer to value function. This will be used to make the old value function V_old
    !accessible to the subroutine used in NAG maximization.
    real(dp), dimension(:,:,:,:,:), pointer :: V_old_pntr

    !optimal index of what_binds (this is used in Howard acceleration algorithm)
    integer, intent(inout) :: optim_wb_gp

    !Values of state variables as saved in the grid (hence the suffix st).
    real(dp), dimension(I_par) :: a_st
    real(dp), dimension(I_par-1) :: rho_st
    real(dp) :: t_st,g_st
    integer :: g_st_ind
    !prepared bounds...
    real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat
    real(dp), dimension(M_par,1) :: rho_min_vec,rho_max_vec

    real(dp), dimension(1,I_par) :: theta_lvl
    real(dp), dimension(M_par,I_par) ::  g_lvl_mat,max_cons_mat
    real(dp), dimension(M_par,M_par) :: P_onerow

    !guess_ind used to cycle over elements of guess_all_gp
    integer :: guess_ind,best_guess_ind

    integer :: m_ind !for cycling over states

    !Variables used in NAG subroutine E04JYF
    real(dp), dimension(ncvar_par) :: xc !choice variables (M_par x I_par matrix in row form)
    real(dp) :: fc !function value in maximization
    !The following two arrays are arrays of reals and integers which can be used to supply information
    !to the maximized function as an alternative to COMMON block.
    real(dp), dimension(1) :: ruser !Contains no information, only defined for purposes of using a NAG subroutine
    integer, dimension(3) :: iuser !will be used to pass what_binds to some subroutines.
    integer :: ifail
    !Bounds:
    real(dp), dimension(M_par*I_par-1) :: bl = 0.0001_dp,bu
    real(dp) :: fmax
    !some other things needed by  e04jyf
    !Try changing the dimension of these things and see what happens (there is only a lower bound)
    integer :: LIW = M_par*I_par +1,LW=max(ncvar_par*(ncvar_par-1)/2+12*ncvar_par,13)
    integer, dimension(M_par*I_par +1) :: IW
    real(dp), dimension(max(ncvar_par*(ncvar_par-1)/2+12*ncvar_par,13)) :: W

    !TMP variables needed for testing E04jcf (remove/comment if I don't use that subroutine
    !in the end)
    integer :: nf
    !Local
    logical :: c_im_found
    real(dp) :: IG_value !Value of initial guess (current + continuation)

    !logical :: debug_tmp_var

    integer :: optim_subroutine = 1 !Temporarily place here,
    !later on put into parameter file. 1 means e04jcf (works on HPC),
    !2 means e04jyf


    !variables used in global optimization
    real(dp), dimension(200) :: comm
    integer :: lcomm = 200 !dimension of comm
    real(dp), dimension(M_par*I_par-1,10) :: list !Second dimension will have to be increased if iinit = 2 or 3
    integer :: sdlist = 10 !This will have to be increased if iinit = 2 or 3
    integer, dimension(M_par*I_par-1) :: numpts = 3 !not actually used in a default initialization method
    integer, dimension(M_par*I_par-1) :: initpt = 5 !not actually used in a default initialization method

    !Variables used in 1-D optimization
    real(dp) :: E1,E2,A,B
    integer :: maxcal

    !loss - used for reporting statistics on loss resulting from violating constraints
    real(dp) :: loss !Loss reported by the maximization program
    real(dp), intent(out) :: loss_gp !loss associated with the optimal choice at the grid point.

    real(dp) :: tmp_val

    !real(dp), dimension(M_par,I_par) :: c_tmp_debug

!    !DEBUG: the following things serve to ascertain difference of outcome when
!    !global optimization used and when local optimization used.
!    real(dp), dimension(ncvar_par) :: xc2 !copy of initial guess
!    real(dp) :: fmax2 !copy of max value

    !(IMP): when calling this, it should be checked that all pointers used in this function
    !and passed using common block are associated - otherwise we will end up with a segmentation fault

    !Common block to be accessible by subroutine CPC_ret
    common /cpc_ret_bloc/ a_st,theta_lvl,rho_st,t_st,g_st,g_st_ind
    !(EFE) : The bounds could actually be prepared even before this subroutine is called,
    !they are the same at every gridpoint.
    common /cpc_ret_bloc2/ g_lvl_mat, P_onerow,l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec

    common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
    !to some subroutines.

    common /commloss/ loss

    optim_subroutine = pars%opt_subroutine

    !If the number of free variables in the optimization problem is 1, then
    !we can't use subroutine e04jcf. This will be the case only when we are solving
    !deterministic version of the model with two players.
    if(M_par*I_par-1 == 1) then
        optim_subroutine = 3
    end if

    !Much of what follows is similar to what was done elsewhere (LFFC).
    !It is advisable to change some of these bits to subroutines or functions so that
    !if we change something, it will be changed everywhere in the program.

    !Get current state from grid:
    a_st(1) = grids%a1_gr(ind_unfold(1))
    a_st(2) = grids%a2_gr(ind_unfold(2))
    rho_st = grids%rho_gr(ind_unfold(3))
    t_st = grids%t_gr(ind_unfold(4))
    g_st = pars%S(1,ind_unfold(5))
    g_st_ind = ind_unfold(5)

    !prepare bounds for checking constraints
    call prep_bounds(pars%l_max,pars%b_min,pars%b_max,pars%a_min,pars%a_max,pars%rho_min,pars%rho_max,&
    l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec)

!        !debug - display bounds
!        if(this_image()==1) then
!        write(*,*) 'bounds before CPC_ret called'
!        write(*,*) 'b_min_mat = ',b_min_mat
!        write(*,*) 'b_max_mat = ',b_max_mat
!        write(*,*) 'l_max_mat = ',l_max_mat
!        write(*,*) 'a_min_mat = ',a_min_mat
!        write(*,*) 'a_max_mat = ',a_max_mat
!        write(*,*) '________________________________'
!        end if

    !Get productivity in levels (the same for all what_binds and consumption choices)
    call prod_lvl(t_st,pars%theta_0,pars%xi_theta,pars%static_environment,theta_lvl)

    !Need government expenditure in levels for all shock realizations. Actually save it in a MxI matrix
    !where every row contains gov. exp. in level in state m
    do m_ind = 1,M_par
        g_lvl_mat(m_ind,:) = pars%S(1,m_ind)*trend_g(pars%mass,theta_lvl,pars%l_max,pars%k_gov)
    end do

    !P_onerow = matmul(reshape([1.0_dp,1.0_dp],[M_par,1]),reshape(pars%P(g_st_ind,:),[1,M_par]))
    P_onerow = pars%P_onerow(g_st_ind,:,:)

    !If we are skipping maximization at this point, compute the 'optimal' value using previous choice
    !and return (Howard).
    if(skip_max) then
        !iuser saves a copy of what_binds -> it is then passed to the
        !NAg subroutines directly (as opposed to a common block).
        iuser = c_guess_all_gp%c_guess_all(optim_wb_gp)%what_binds

        !Transform the initial guess into vector form so that it can be used
        !by NAG subroutines.
        call c_inc_mtv(c_guess_all_gp%c_guess_all(optim_wb_gp)%c_guess,xc)

        !Compute the value of initial guess.
        call CPC_ret(ncvar_par,xc,IG_value,iuser,ruser)
        !Save the value of the initial guess as the new value function. Initial guess is not
        !updated (i.e., policy function the same)
        V_new_gp = -IG_value

        return
    end if

    call max_feas_cons(pars%l_max,pars%mass,theta_lvl,g_lvl_mat,max_cons_mat)
    !Use this to generate the vector of upper bounds...

    !get upper bounds
    call c_inc_mtv(max_cons_mat,bu)

    V_new_gp = -10.0_dp**20 !initialize value (must be very low to make sure that this
    !does not serve as an artificial lower bound and also so we detect errors
    !arising from its usage easily). Also if at some point we didn't
    !find a value better than this, then the optimal policy would not be updated.


    if(glob_opt) then !use global optimization
        !Call initialization routing (has to be done only once)
        call e05jaf(0,comm,lcomm,ifail)
    end if

    !Need to cycle over elements of what_binds. Get a maximum at every index.
    !If there is -1, skip to the next row of what_binds...
    !update c_guess for every what_binds and save optimal value. Return the
    !highest value of those found.
    do guess_ind = 1,size(c_guess_all_gp%c_guess_all)
        !If the first element of what_binds is -1, it means that for this particular
        !element of c_guess_all_gp, no feasible solution was found.
        !Go to the next guess.
        if(c_guess_all_gp%c_guess_all(guess_ind)%what_binds(1) == -1) go to 100

        !debug_tmp_var = .false. !so this will be true if we skipped the cycle every time...

        !iuser saves a copy of what_binds -> it is then passed to the
        !NAg subroutines directly (as opposed to a common block).
        iuser = c_guess_all_gp%c_guess_all(guess_ind)%what_binds

        !Transform the initial guess into vector form so that it can be used
        !by NAG subroutines.
        call c_inc_mtv(c_guess_all_gp%c_guess_all(guess_ind)%c_guess,xc)
        !need to pass what_binds -> as iuser? Or as common...?

        !Compute the value of initial guess. This will be a lower bound in case
        !the maximization subroutine returns something crazy
        call CPC_ret(ncvar_par,xc,IG_value,iuser,ruser)

        loss_gp = loss !save the loss at the initial guess (unless something better found later,
        !this will be the point chosen and the associated loss will be the loss at this gp)

        ifail = 1 !silent exit
        !Find maximum. In the first few iterations, global optimization is used (due to
        !non-convexity of feasible sets). In later iterations, once we have a feasible
        !guess hopefully in the convex subset of the feasible set in which the global optimum is,
        !we use local optimization subroutines - these are much faster and also
        !usually give better results for good initial guess.

        if(glob_opt) then !use global optimization
              !experiment with various iinit values (the input to the left of bl)
              !xc is not actually used as an initial guess if we used
              !global optimization - this could be somehow used only if
              !iinit is > 1? So can save a LOT of time by not getting initial guess
              !in the first iteration (if global maximization is used).
              call e05jbf(M_par*I_par-1,CPC_ret3,0,0,bl,bu,sdlist,list,numpts,initpt,E05JBK, &
                xc,fmax,comm,lcomm,iuser,ruser,ifail)
                !if ifail>0, then fmax could contain some nonsense (and xc could be nonsense).
                !If this happens, perhaps use one of the local minimization subroutines
                !starting at the initial guess...
        else

        select case(optim_subroutine)
            case(1) !Quadratic approximation method
                call e04jcf(CPC_ret2,M_par*I_par-1,ncvar_par*2+1,xc,bl,bu,0.1_dp,&
                pars%rhoend,e04jcp,pars%VFI_opt_maxiter,fmax,nf,iuser,ruser,ifail)
            case(2) !Newton-like method
                call e04jyf(M_par*I_par-1,0,CPC_ret,bl,bu,xc,fmax,iw,liw,w,lw,iuser,ruser,ifail)
            case(3) !1-D optimization (only used then M_par = 1, I_par = 2)
                A = bl(1)
                B = bu(1)
                e1 = pars%e1
                e2 = pars%e2
                maxcal = pars%VFI_opt_maxiter

                call E04ABA(CPC_ret_1d,E1,E2,A,B,maxcal,XC(1),fmax,IUSER,RUSER,IFAIL)

            case default
                write(*,*) 'Error in subroutine VFI_solve_gp - wrong value of optim_subroutine'
                error stop
        end select


        end if


        !Check whether the optimum found by one of the subroutines is better than the
        !value of the initial guess. If yes, update the initial guess and use its
        !value as the new value. Otherwise do not update the initial guess and
        !use its value.
        if(fmax<IG_value) then !(optimum found is better than init. guess, remember we work with (-1) multiples)
            !Update initial guess (have to reshape it into incomplete cons. matrix first)
            call c_inc_vtm(xc,c_guess_all_gp%c_guess_all(guess_ind)%c_guess)
            !If the maximum found here is better than all other maxima
            !found at this gridpoint (for other values of what_binds),
            !save the optimal value and choice into value function
            !and policy function (can still be rewritten if better things found for
            !different what_binds).
            if(-fmax>V_new_gp) then
                V_new_gp = -fmax
                C_pol_gp=c_guess_all_gp%c_guess_all(guess_ind)%c_guess !the right side already contains updated guess
                !Need to get the last elemement separately.
                c_im_found = .false.
                call get_c_im(C_pol_gp,c_guess_all_gp%c_guess_all(guess_ind)%what_binds,&
                a_st,rho_st,theta_lvl,g_lvl_mat(:,1),g_st_ind,pars,P_onerow,c_im_found)

                if(.not. c_im_found) then
                    call fix_c(C_pol_gp,rho_st,pars)
                    !write(*,*) 'using fix_c in subroutine VFI_solve_gp (place 1)'
                    !This should be very rare (only in cases where the initial guess was absolutely terrible)
                    !Recompute the value and save it (we assume that fixing the consumption from being negative
                    !improved the value - this should be the case because value is terribly penalized even for tiny
                    !violations of the non-negativity constraint on consumption).
                    call c_inc_mtv(C_pol_gp,xc)
                    call CPC_ret(ncvar_par,xc,tmp_val,iuser,ruser)
                    V_new_gp = -tmp_val
                    c_guess_all_gp%c_guess_all(guess_ind)%c_guess = C_pol_gp

                end if

                !Also save the index of the optimal choice:
                optim_wb_gp = guess_ind

                if(glob_opt) then
                    !The last time a global optimization subroutine called a CPC_ret subroutine
                    !does not necessarily correspond to the optimum - so for the optimum choice we have to
                    !(if we care enough) compute the loss again.
                    loss_gp = 0.0_dp
                else
                    loss_gp = loss !save loss (this is passed using a common block
                    !from subroutine CPC_ret)
                end if
            end if
        else !Use initial guess instead (and don't update initial guess)
            if(-IG_value>V_new_gp) then
                V_new_gp = -IG_value
                C_pol_gp=c_guess_all_gp%c_guess_all(guess_ind)%c_guess !here the right side not updated
                !Need to get the last elemement separately.
                c_im_found = .false.
                call get_c_im(C_pol_gp,c_guess_all_gp%c_guess_all(guess_ind)%what_binds,&
                a_st,rho_st,theta_lvl,g_lvl_mat(:,1),g_st_ind,pars,P_onerow,c_im_found)
                if(.not. c_im_found) then
                    call fix_c(C_pol_gp,rho_st,pars)
                    write(*,*) 'using fix_c in subroutine VFI_solve_gp (place 2)'
                    !In this case the initial guess was better than what we found even if one of its
                    !components had negative consumption. Fixing this only improves the value, therefore
                    !we do not need to check whether this point should be used (because even slightly
                    !negative consumption was hugely penalized). Nevertheless, we update the value
                    !and the fixed guess.
                    call c_inc_mtv(C_pol_gp,xc)
                    call CPC_ret(ncvar_par,xc,tmp_val,iuser,ruser)
                    V_new_gp = -tmp_val
                    c_guess_all_gp%c_guess_all(guess_ind)%c_guess = C_pol_gp
                end if
                !Save optimal wb index (optimal in the sense that even though maximization
                !failed to find something good, initial guess for this grid point was the
                !best).
                optim_wb_gp = guess_ind
            end if

        end if

100 end do

!if(V_new_gp < -10.0_dp**10) then
!    write(*,*) 'Something bad happened...'
!    error stop
!end if
end subroutine VFI_solve_gp




!function theta_prime computes next-period state. Input is current productivity
!of all agents in levels.
function theta_prime(theta_lvl,theta_0,xi_theta,static_environment,theta_gr)
    real(dp) :: theta_prime
    real(dp), dimension(1,I_par), intent(in) :: theta_lvl, theta_0, xi_theta
    logical, intent(in) :: static_environment
    real(dp), dimension(:), intent(in) :: theta_gr

    real(dp), dimension(1,I_par) :: theta_prime_lvl


    if(static_environment) then
        theta_prime_lvl = theta_0 !could just take theta_lvl, the problem is
        !that there could be compounded numerical errors over many iterations
        !so it's safer to do this
    else
        theta_prime_lvl = theta_lvl * (xi_theta + 1.0_dp)
    end if

    theta_prime = theta_prime_lvl(1,2)/theta_prime_lvl(1,1)

    !If we got outside of the grid, truncate. This should happen after
    !many periods as long as the grid is sufficiently large.
    theta_prime = max(theta_prime,minval(theta_gr))
    theta_prime = min(theta_prime,maxval(theta_gr))

end function theta_prime

!Subroutine solve_FP_problem solves the first period problem (t=0), given a value function
!for periods t>=1 and initial conditions.
!
!It is better to have one subroutine for solving the FP problem and one for solving
!the subsequent periods because the dimensionality of choice variables is different
!and there are a lot of other small changes which would make it confusing if we
!put all simulation into one subroutine).
subroutine solve_FP_problem(V_old,b_init,theta_0,g_ind_init,pars,FP_result,t_gr)
    Use nag_library, Only: e04jyf, nag_wp, e04jcf, e04jcp,e05jaf,e05jbf,E05JBK

    real(dp), dimension(:,:,:,:,:), intent(in) :: V_old
    real(dp), dimension(1,I_par), intent(in) :: b_init
    real(dp), dimension(1,I_par), intent(in) :: theta_0 !initial productivity (in levels)
    integer, intent(in) :: g_ind_init
    type(par), intent(in) :: pars
    type(FP_result_type), intent(out) :: FP_result
    real(dp), dimension(:), intent(in) :: t_gr

    !(For now just use save_series = .true. Later on can return them as arguments,
    !this will be useful if we want to simulate many series and compare results, etc.)

    !Some of these variables may not be needed in the end, check this at the end and delete the
    !unnecessary ones.
    real(dp), dimension(M_par,I_par) ::  g_lvl_mat,max_cons_mat!, P_onerow

    !Variables used in NAG subroutine E04JYF
    real(dp), dimension(I_par) :: xc,xc2,xc_tmp !choice variables are cons of players 1,...,I
    real(dp) :: fc !function value in maximization
    !The following two arrays are arrays of reals and integers which can be used to supply information
    !to the maximized function as an alternative to COMMON block.
    real(dp), dimension(1) :: ruser !Contains no information, only defined for purposes of using a NAG subroutine
    integer, dimension(3) :: iuser !will be used to pass what_binds to some subroutines.
    integer :: ifail
    !Bounds:
    real(dp), dimension(I_par) :: bl = 0.0001_dp,bu
    real(dp) :: fmax,fmax2,fmax_tmp
    !some other things needed by  e04jyf
    !Try changing the dimension of these things and see what happens (there is only a lower bound)
    integer :: LIW = I_par +2,LW=max(I_par*(I_par-1)/2+12*I_par,13)
    integer, dimension(I_par + 2) :: IW
    real(dp), dimension(max(I_par*(I_par-1)/2+12*I_par,13)) :: W

    integer :: max_ind

    !TMP variables needed for testing E04jcf (remove/comment if I don't use that subroutine
    !in the end)
    integer :: nf
    !Local
    logical :: c_im_found
    real(dp) :: IG_value !Value of initial guess (current + continuation)

    !prepared bounds...
    real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat
    real(dp), dimension(M_par,1) :: rho_min_vec,rho_max_vec

    !variables used in global optimization
    real(dp), dimension(200) :: comm
    integer :: lcomm = 200 !dimension of comm
    real(dp), dimension(I_par,10) :: list !Second dimension will have to be increased if iinit = 2 or 3
    integer :: sdlist = 10 !This will have to be increased if iinit = 2 or 3
    integer, dimension(I_par) :: numpts = 3 !not actually used in a default initialization method
    integer, dimension(I_par) :: initpt = 5 !not actually used in a default initialization method

    integer :: m_ind
    integer :: tmp_test_var !remove this later

    !Copies of data which we will pass to CPC subroutine using common block (need
    !a copy because these object are passed as arguments to the subroutine
    !so they can't be in a common block).
    real(dp), dimension(1,I_par) :: b_init_cp
    real(dp), dimension(1,I_par) :: theta_0_cp
    integer :: g_ind_init_cp

        !The following variables are use for randomly drawing initial shock realization
        !if s_init_index = -1 in parameter file. We want to draw the shock with frequency associated
        !with stationary distribution.
        real(dp) :: k,r1,r2,randnum


    !Debug - possibly remove later
    real(dp), dimension(1,I_par) :: c_ss,l_ss

    !Debug values
    real(dp) :: fc_debug
    real(dp), dimension(I_par) :: xc_debug
    integer :: inform_debug

    integer :: attmaxloc


    !Declare the common block use in CPC return computation.
    common /cpc_ret_fp/ b_init_cp, theta_0_cp,g_ind_init_cp

    !common /cpc_ret_fp_bnd/ l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec
    common /cpc_ret_bloc_fp/ g_lvl_mat,l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec

    !Create copies of inputs which need to be made available to subroutine
    !called by the maximization subroutine using a common block.
    b_init_cp = b_init
    theta_0_cp = theta_0
    g_ind_init_cp = g_ind_init

    !Save the actual index in FP_result!
    FP_result%s_ind_init = g_ind_init_cp

    !Need government expenditure in levels for all shock realizations. Actually save it in a MxI matrix
    !where every row contains gov. exp. in level in state m. Here we really need only one row
    do m_ind = 1,M_par
        g_lvl_mat(m_ind,:) = pars%S(1,m_ind)*trend_g(pars%mass,theta_0_cp,pars%l_max,pars%k_gov)
    end do

    !The following - should be moved to a higher level - they are already generated when VFI_solve_GP
    !is called for the last time...
    !prepare bounds for checking constraints
    call prep_bounds(pars%l_max,pars%b_min,pars%b_max,pars%a_min,pars%a_max,pars%rho_min,pars%rho_max,&
    l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec)

    !Generate bounds ul and ub
    bl = 0.000001_dp

    call max_feas_cons(pars%l_max,pars%mass,theta_0_cp,g_lvl_mat,max_cons_mat)
    bu = max_cons_mat(g_ind_init_cp,:) !depends on the current shock realization


    !DEBUG: check inputs of cpc_ret2 and ret3


    !Initialization routine for e05jbf
    call e05jaf(0,comm,lcomm,ifail)

    !Need to have bounds prepared...... And call the initialization subroutine...
    ifail = -1 !set 1 for silent exit
    numpts = 100
    !call e05jbf(I_par,CPC_ret3_FP,0,0,bl,bu,sdlist,list,numpts,initpt,E05JBK, &
    !xc,fmax,comm,lcomm,iuser,ruser,ifail)

    !Run global optimization subroutine with lots of different initial conditions
    !(they are random) - we can afford this robustness in the first period
    fmax_tmp = 100000000000000 !initialize (working with minus 1 multiples so lower number is better here)
    xc_tmp = 0.0000000001_dp

    !If perm_eps = 0, we can afford to try more times.
    if(pars%perm_eps == 0.0_dp) then
        attmaxloc = pars%sim_attmax1 * 5
    else
        attmaxloc = pars%sim_attmax1
    end if

    do max_ind = 1,attmaxloc
        ifail = 1
        call e05jbf(I_par,CPC_ret3_FP,0,4,bl,bu,sdlist,list,numpts,initpt,E05JBK, &
            xc,fmax,comm,lcomm,iuser,ruser,ifail)
        if(fmax<fmax_tmp) then
            fmax_tmp = fmax
            xc_tmp = xc
        end if
    end do
    xc = xc_tmp
    fmax = fmax_tmp

    !Now use the optimum choice found by the globalization subroutine to
    !initialize search for a local optimum (could be an improvement).
    xc2 = xc
    !test: try to use steady state as initial guess here

    call e04jcf(CPC_ret2_FP,I_par,I_par*2+1,xc2,bl,bu,0.1_dp,&
    0.000001_dp,e04jcp,5000,fmax2,nf,iuser,ruser,ifail)

    !This thing is largely redundant now that we use global optimization many times with
    !random guesses - but it indicate instability and reveal some bugs so I'm gonna
    !leave it here...

    if(abs(fmax2-fmax) > 0.0001_dp) then
        write(*,*) 'Warning: The local and the global maximization in first period yield different results.'
        write(*,*) 'Difference in function value = ', (fmax2-fmax)
    end if

    !Take the best of these results (it is the one with the lower value since
    !we minimized a (-1) multiple)
    if(fmax<fmax2) then
        FP_result%c = reshape(xc,[1,I_par])
    else
        FP_result%c = reshape(xc2,[1,I_par])
    end if


    !Compute the remainder of FP_result like in subroutine CPC_ret3_FP (I should write a subroutine for
    !this at some point)
    !Given consumption, compute labour supply (using resource constraint and intratemporal EE)
    call lab_FP(FP_result%c,FP_result%l,theta_0_cp,pars%mass,g_lvl_mat(g_ind_init_cp,1),pars%Gamma_par)

    !Get next-period asset holdings:
    call b_prime_fp(b_init_cp,FP_result%c,FP_result%l,pars%A_par,pars%B_par,pars%Gamma_par,FP_result%b_prime)

    !Get marginal utility of consumption
    call util_c(FP_result%c,FP_result%u_c,pars%A_par)
    !get a_prime
    FP_result%a_prime = (FP_result%b_prime*FP_result%u_c)/pars%beta

    !rho_prime: This will need to be generalized for I>2
    FP_result%rho_prime = FP_result%u_c(1,1)/FP_result%u_c(1,2)

    !Get next-period 'time' (assumes that the first point of the grid is the time in the
    !second period (index t=1)
    FP_result%t_prime = t_gr(1)



    !get marginal utility of labour
    call util_l(FP_result%l,FP_result%u_l,pars%B_par,pars%Gamma_par)



end subroutine solve_FP_problem

!Subroutine CPC_ret3_FP computes the current + continuation return
!in first-period (t=0). It is written for use with a global optimization subroutine.
subroutine CPC_ret3_FP(n,xc,fc,nstate,iuser,ruser,inform)
        real (dp), intent (out) :: fc
        integer, intent (in) :: n
        integer, intent(in) :: nstate !=1 if the subroutine is called for the first time,
        !not actually used anywhere atm but required by NAg.
        real (dp), intent (inout) :: ruser(*)
        real (dp), intent (in) :: xc(n)
        integer, intent (inout) :: iuser(*) !In this instance, it is used to pass what_binds
        integer, intent(out) :: inform

        real(dp), dimension(:,:,:,:,:), pointer :: V_old_pntr !pointer to V_old, passed using common block

        !The following variables are passed to the function using a common block
        real(dp), dimension(:), pointer :: a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr
        real(dp), dimension(:,:), pointer :: S_pntr,P_pntr, mass_pntr
        !and transition matrix.
        real(dp), pointer :: A_par_pntr,B_par_pntr,Gamma_par_pntr
        type(par), pointer :: pars_pntr

        !Values of state variables as saved in the grid (hence the suffix st).
        real(dp), dimension(I_par) :: a_st
        real(dp), dimension(I_par-1) :: rho_st
        real(dp) :: t_st,g_st
        integer :: g_st_ind
        real(dp), dimension(1,I_par) :: theta_lvl

        real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat
        real(dp), dimension(M_par,1) :: rho_min_vec,rho_max_vec

        real(dp), dimension(M_par,I_par) :: g_lvl_mat
        real(dp), dimension(M_par,M_par) :: P_onerow

        !Local variables
        real(dp), dimension(1,I_par) :: c_vec, u_c_vec = 0.0_dp !consumption vector and MU_c vector
        real(dp), dimension(1,I_par) :: U !utility of both agents
        !Labour
        real(dp), dimension(1,I_par) :: lab_vec = 0.0_dp !labour vector
        !next-period states
        real(dp) :: t_st_pr
        real(dp), dimension(1,I_par) :: a_prime,b_prime
        real(dp), dimension(1,1) :: rho_prime !with more states it's a column vector

        logical :: c_im_found = .false., constr_ok

        integer :: g_ind !index of shock realization

        real(dp) :: loss !loss from violating constraints.

        real(dp) :: conval !continuation value (pre-discount)

        real(dp), dimension(M_par,1) ::  ones = 1.0_dp

        !Current state will be passed as a common block
        real(dp), dimension(1,I_par) :: b_init_cp
        real(dp), dimension(1,I_par) :: theta_0_cp
        integer :: g_ind_init_cp

        real(dp) :: trend_g

        common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
        !to some subroutines.
        common /pars_pointer/ pars_pntr

        common /cpc_ret_bloc/ a_st,theta_lvl,rho_st,t_st,g_st,g_st_ind
        common /cpc_ret_bloc_fp/ g_lvl_mat,l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec

        !Need to pass state to this thing.

        common /comm1/ a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr,A_par_pntr,B_par_pntr,&
        Gamma_par_pntr,S_pntr,P_pntr,mass_pntr

        common /cpc_ret_fp/ b_init_cp, theta_0_cp,g_ind_init_cp

        !Negative value lead to termination of maximization procedure.
        inform = 0

        c_vec = reshape(xc,[1,I_par]) !Convenient to have this so we don't have to use
        !reshape all the time (element (1,i) contains consumption of agent i in period 0)

        !Given consumption, compute labour supply (using resource constraint and intratemporal EE)
        call lab_FP(c_vec,lab_vec,theta_0_cp,mass_pntr,g_lvl_mat(g_ind_init_cp,1),Gamma_par_pntr)

        !Get next-period asset holdings:
        call b_prime_fp(b_init_cp,c_vec,lab_vec,A_par_pntr,B_par_pntr,Gamma_par_pntr,b_prime)

        !Get marginal utility of consumption
        call util_c(c_vec,u_c_vec,A_par_pntr)

        !get a_prime
        a_prime = (b_prime*u_c_vec)/pars_pntr%beta

        !This will need to be generalized for I>2
        rho_prime(1,1) = u_c_vec(1,1)/u_c_vec(1,2)

        !Get t_st_pr (next-period 'time')
        t_st_pr = t_gr_pntr(1)

        call util(c_vec,lab_vec,U,A_par_pntr,B_par_pntr,Gamma_par_pntr) !need to use mod_util for this



        !Now we have all current choice variables and next-period states -> we can
        !compute current + continuation return

        !Current return
        fc = sum(pars_pntr%mass*pars_pntr%alpha*U)


        !Continuation return
        !1) Get the loss associated with violating constraints
        !!!The bounds prepared here were in fact prepared earlier by VFI_solve_GP
        !and end up here because they're in a common block.

        !Here we're got state-independent things - stack the rows in a matrix M_par times so
        !we can use the subroutine for checking constraints on state-dependent plans,
        !then divide by M_par. It's less cumbersome than writing a special subroutine
        !or generalizing check_constr. Maybe change later, not a priority atm.
        call check_constr(matmul(ones,c_vec),matmul(ones,lab_vec),matmul(ones,b_prime),&
        b_min_mat,b_max_mat,l_max_mat,theta_0_cp,matmul(ones,a_prime),a_min_mat&
        ,a_max_mat,iuser,rho_min_vec,rho_max_vec,constr_ok,loss,pars_pntr%A_par,.false.,.false.)

        !Don't forget to divide loss by M_par (so it's not counted here more than in later periods)
        fc = fc - loss/real(M_par,dp)


        !2) Truncate (instead of extrapolation, use nearest neighbour extrap + loss)
        !(just like in CPC_ret and other functions but have to use reshape here
        !because the constraints are in matrix form
        if(.not. constr_ok) then
            a_prime=max(reshape(a_min_mat(1,:),[1,2]),a_prime)
            a_prime=min(reshape(a_max_mat(1,:),[1,2]),a_prime)
            rho_prime=max(reshape(rho_min_vec(1,:),[1,1]),rho_prime)
            rho_prime=min(reshape(rho_max_vec(1,:),[1,1]),rho_prime)
        end if

        !3) Use interpolation to get the continuation value
        call interp_V(V_old_pntr,[a_prime(1,1),a_prime(1,2),rho_prime(1,1),t_st_pr]&
        ,g_ind_init_cp,conval,a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr,size(a1_gr_pntr),&
        size(rho_gr_pntr),size(t_gr_pntr))

        fc = fc + pars_pntr%beta*conval!+ discounted continuation value (loss vas already subtracted and is not discounted)

        !At the end multiply fc by -1 (because we want to use it in a minimization subroutine)
        fc = fc*(-1.0_dp)
end subroutine CPC_ret3_FP


!Subroutine CPC_ret2_FP is a copy of CPC_ret3 with a different interface so it can
!be used by a different NAg subroutine - really should just write different interface
!subroutines which would call a single subroutine which would do all the computation
!Helpful for clarity of code and avoiding bugs.
subroutine CPC_ret2_FP(n,xc,fc,iuser,ruser,inform)
        real (dp), intent (out) :: fc
        integer, intent (in) :: n
        real (dp), intent (inout) :: ruser(*)
        real (dp), intent (in) :: xc(n)
        integer, intent (inout) :: iuser(*) !In this instance, it is used to pass what_binds
        integer, intent(out) :: inform

        real(dp), dimension(:,:,:,:,:), pointer :: V_old_pntr !pointer to V_old, passed using common block

        !The following variables are passed to the function using a common block
        real(dp), dimension(:), pointer :: a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr
        real(dp), dimension(:,:), pointer :: S_pntr,P_pntr, mass_pntr
        !and transition matrix.
        real(dp), pointer :: A_par_pntr,B_par_pntr,Gamma_par_pntr
        type(par), pointer :: pars_pntr

        !Values of state variables as saved in the grid (hence the suffix st).
        real(dp), dimension(I_par) :: a_st
        real(dp), dimension(I_par-1) :: rho_st
        real(dp) :: t_st,g_st
        integer :: g_st_ind
        real(dp), dimension(1,I_par) :: theta_lvl

        real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat
        real(dp), dimension(M_par,1) :: rho_min_vec,rho_max_vec

        real(dp), dimension(M_par,I_par) :: g_lvl_mat
        real(dp), dimension(M_par,M_par) :: P_onerow

        !Local variables
        real(dp), dimension(1,I_par) :: c_vec, u_c_vec = 0.0_dp !consumption vector and MU_c vector
        real(dp), dimension(1,I_par) :: U !utility of both agents
        !Labour
        real(dp), dimension(1,I_par) :: lab_vec = 0.0_dp !labour vector
        !next-period states
        real(dp) :: t_st_pr
        real(dp), dimension(1,I_par) :: a_prime,b_prime
        real(dp), dimension(1,1) :: rho_prime !with more states it's a column vector

        logical :: c_im_found = .false., constr_ok

        integer :: g_ind !index of shock realization

        real(dp) :: loss !loss from violating constraints.

        real(dp) :: conval !continuation value (pre-discount)

        real(dp), dimension(M_par,1) ::  ones = 1.0_dp

        !Current state will be passed as a common block
        real(dp), dimension(1,I_par) :: b_init_cp
        real(dp), dimension(1,I_par) :: theta_0_cp
        integer :: g_ind_init_cp

        real(dp) :: trend_g




        common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
        !to some subroutines.
        common /pars_pointer/ pars_pntr

        common /cpc_ret_bloc/ a_st,theta_lvl,rho_st,t_st,g_st,g_st_ind
        common /cpc_ret_bloc_fp/ g_lvl_mat,l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec

        !Need to pass state to this thing.

        common /comm1/ a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr,A_par_pntr,B_par_pntr,&
        Gamma_par_pntr,S_pntr,P_pntr,mass_pntr

        common /cpc_ret_fp/ b_init_cp, theta_0_cp,g_ind_init_cp

        !Negative value lead to termination of maximization procedure.
        inform = 0

        c_vec = reshape(xc,[1,I_par]) !Convenient to have this so we don't have to use
        !reshape all the time (element (1,i) contains consumption of agent i in period 0)

        !DEBUG!!! remove later or the program will not work!!!
        !c_vec = 0.286730230002313_dp !the steady-state value

        !Given consumption, compute labour supply (using resource constraint and intratemporal EE)
        call lab_FP(c_vec,lab_vec,theta_0_cp,mass_pntr,g_lvl_mat(g_ind_init_cp,1),Gamma_par_pntr)


        !Get next-period asset holdings:
        call b_prime_fp(b_init_cp,c_vec,lab_vec,A_par_pntr,B_par_pntr,Gamma_par_pntr,b_prime)

        !Get marginal utility of consumption
        call util_c(c_vec,u_c_vec,A_par_pntr)

        !get a_prime
        a_prime = (b_prime*u_c_vec)/pars_pntr%beta

        !This will need to be generalized for I>2
        rho_prime(1,1) = u_c_vec(1,1)/u_c_vec(1,2)

        !Get next-period 'time'
        t_st_pr = t_gr_pntr(1)

        call util(c_vec,lab_vec,U,A_par_pntr,B_par_pntr,Gamma_par_pntr) !need to use mod_util for this



        !Now we have all current choice variables and next-period states -> we can
        !compute current + continuation return

        !Current return
        fc = sum(pars_pntr%mass*pars_pntr%alpha*U)


        !Continuation return
        !1) Get the loss associated with violating constraints
        !!!The bounds prepared here were in fact prepared earlier by VFI_solve_GP
        !and end up here because they're in a common block.

        !Here we're got state-independent things - stack the rows in a matrix M_par times so
        !we can use the subroutine for checking constraints on state-dependent plans,
        !then divide by M_par. It's less cumbersome than writing a special subroutine
        !or generalizing check_constr. Maybe change later, not a priority atm.
        call check_constr(matmul(ones,c_vec),matmul(ones,lab_vec),matmul(ones,b_prime),&
        b_min_mat,b_max_mat,l_max_mat,theta_0_cp,matmul(ones,a_prime),a_min_mat&
        ,a_max_mat,iuser,rho_min_vec,rho_max_vec,constr_ok,loss,pars_pntr%A_par,.false.,.false.)

        !Don't forget to divide loss by M_par (so it's not counted here more than in later periods)
        fc = fc - loss/real(M_par,dp)

        !2) Truncate (instead of extrapolation, use nearest neighbour extrap + loss)
        !(just like in CPC_ret and other functions but have to use reshape here
        !because the constraints are in matrix form
        if(.not. constr_ok) then
            a_prime=max(reshape(a_min_mat(1,:),[1,2]),a_prime)
            a_prime=min(reshape(a_max_mat(1,:),[1,2]),a_prime)
            rho_prime=max(reshape(rho_min_vec(1,:),[1,1]),rho_prime)
            rho_prime=min(reshape(rho_max_vec(1,:),[1,1]),rho_prime)
        end if


        !3) Use interpolation to get the continuation value
        call interp_V(V_old_pntr,[a_prime(1,1),a_prime(1,2),rho_prime(1,1),t_st_pr]&
        ,g_ind_init_cp,conval,a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr,size(a1_gr_pntr),&
        size(rho_gr_pntr),size(t_gr_pntr))

        !DEBUG
        !write(*,*) 'conval = ',conval
        !write(*,*) 'disc. conval = ',pars_pntr%beta*conval

        fc = fc + pars_pntr%beta*conval!+ discounted continuation value (loss vas already subtracted and is not discounted)

        !At the end multiply fc by -1 (because we want to use it in a minimization subroutine)
        fc = fc*(-1.0_dp)

end subroutine CPC_ret2_FP

!Subroutine sim_result_alloc allocates memory to allocatable arrays of
!sim_result_type. It contains no checks of whether the arrays are already allocated
!so it needs to be used carefully, otherwise it will lead to crashes, etc.
subroutine sim_result_alloc(sim_result,pars,wtd)
    type(sim_result_type), intent(inout) :: sim_result
    type(par), intent(in) :: pars
    integer, intent(in) :: wtd !what to do (1 = allocate, 0 = deallocate)

    if(wtd == 1) then
        !Allocate memory to simulated series
    allocate(sim_result%c_sim(pars%T_sim,I_par),sim_result%l_sim(pars%T_sim,I_par),sim_result%b_sim(pars%T_sim,I_par),&
    sim_result%b_prime_sim(pars%T_sim,I_par),sim_result%a_prime_sim(pars%T_sim,I_par),sim_result%u_c_sim(pars%T_sim,I_par),&
    sim_result%u_l_sim(pars%T_sim,I_par),sim_result%theta_lvl_sim(pars%T_sim,I_par))
    allocate(sim_result%rho_prime_sim(pars%T_sim),sim_result%t_prime_sim(pars%T_sim),&
    sim_result%t_sim(pars%T_sim),sim_result%s_ind_sim(pars%T_sim),sim_result%s_real_sim(pars%T_sim),&
    sim_result%R_sim(pars%T_sim-1),sim_result%tau_sim(pars%T_sim),sim_result%g_lvl_sim(pars%T_sim),&
    sim_result%gini_income_sim(pars%T_sim),sim_result%gini_prod_sim(pars%T_sim))
    allocate(sim_result%Y_sim(pars%T_sim),sim_result%B_gov_sim(pars%T_sim),&
    sim_result%B_gov_to_GDP_sim(pars%T_sim))
    allocate(sim_result%G_to_Y_sim(pars%T_sim),sim_result%gini_allincome_sim(pars%T_sim))

    else if(wtd==0) then
        !Deallocate
    deallocate(sim_result%c_sim,sim_result%l_sim,sim_result%b_sim,&
    sim_result%b_prime_sim,sim_result%a_prime_sim,sim_result%u_c_sim,&
    sim_result%u_l_sim,sim_result%theta_lvl_sim)
    deallocate(sim_result%rho_prime_sim,sim_result%t_prime_sim,&
    sim_result%t_sim,sim_result%s_ind_sim,sim_result%s_real_sim,&
    sim_result%R_sim,sim_result%tau_sim,sim_result%g_lvl_sim,&
    sim_result%gini_income_sim,sim_result%gini_prod_sim)
    deallocate(sim_result%Y_sim,sim_result%B_gov_sim,&
    sim_result%B_gov_to_GDP_sim)
    deallocate(sim_result%G_to_Y_sim,sim_result%gini_allincome_sim)
    else
        write(*,*) 'Error: wrong value of wtd in subroutine sim_result_alloc.'
        sync all
        error stop
    end if


end subroutine sim_result_alloc


!Subroutine sim_series computes simulated series for periods 1,...,T_sim-1
!(so that the total number of simulated elements including thhose for t=0
!will be T_sim)
subroutine sim_series(sim_result,FP_result,pars,grids)
    Use nag_library, Only: e04jyf, nag_wp, e04jcf, e04jcp,e05jaf,e05jbf,E05JBK
    !The E05 libraries are not used at this point but may be if we
    !want to look for corner solutions as well (in which case we probably would not have
    !a good initial guess from FP problem and we should use maximization subroutines).

    type(sim_result_type), intent(inout) :: sim_result
    type(FP_result_type), intent(in) :: FP_result
    type(par), intent(in) :: pars
    type(grids_type),intent(in) :: grids

    real(dp), dimension(M_par,I_par) :: c_guess !This will contain guess of consumption
    !that will be passed between periods. So far corner solutions not implemented in simulation
    !(only in computing value function) -> later on change this to passing c_guess_all
    !in between periods, but in that case we will have to peform a search for
    !initial guess in the first period (because we don't have a full guess
    !from the first-period problem - something like LFFC but at one point only).
    real(dp), dimension(M_par,I_par) :: c_guess_old !backup of initial guess for cases where
    !optimization fails for various reasons.

    real(dp), dimension(M_par,1) :: ones_col = 1.0_dp

    integer :: t_ind !period index

    !Pointer to value function. This will be used to make the old value function V_old
    !accessible to the subroutine used in NAG maximization.
    real(dp), dimension(:,:,:,:,:), pointer :: V_old_pntr


    !Values of state variables as saved in the grid (hence the suffix st).
    real(dp), dimension(I_par) :: a_st
    real(dp), dimension(I_par-1) :: rho_st
    real(dp) :: t_st,g_st
    integer :: g_st_ind
    !prepared bounds...
    real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat
    real(dp), dimension(M_par,1) :: rho_min_vec,rho_max_vec

    real(dp), dimension(1,I_par) :: theta_lvl
    real(dp), dimension(M_par,I_par) ::  g_lvl_mat,max_cons_mat
    real(dp) :: t_st_pr
    real(dp), dimension(M_par,M_par) :: P_onerow

    !guess_ind used to cycle over elements of guess_all_gp
    !integer :: guess_ind,best_guess_ind

    integer :: m_ind !for cycling over states

    !Variables used in NAG subroutine E04JYF
    real(dp), dimension(ncvar_par) :: xc !choice variables (M_par x I_par matrix in row form)
    real(dp) :: fc !function value in maximization
    !The following two arrays are arrays of reals and integers which can be used to supply information
    !to the maximized function as an alternative to COMMON block.
    real(dp), dimension(1) :: ruser !Contains no information, only defined for purposes of using a NAG subroutine
    integer, dimension(3) :: iuser !will be used to pass what_binds to some subroutines.
    integer :: ifail
    !Bounds:
    real(dp), dimension(M_par*I_par-1) :: bl = 0.0001_dp,bu
    real(dp) :: fmax
    !some other things needed by  e04jyf
    !Try changing the dimension of these things and see what happens (there is only a lower bound)
    integer :: LIW = M_par*I_par +1,LW=max(ncvar_par*(ncvar_par-1)/2+12*ncvar_par,13)
    integer, dimension(M_par*I_par +1) :: IW
    real(dp), dimension(max(ncvar_par*(ncvar_par-1)/2+12*ncvar_par,13)) :: W

    !TMP variables needed for testing E04jcf (remove/comment if I don't use that subroutine
    !in the end)
    integer :: nf
    !Local
    logical :: c_im_found
    real(dp) :: IG_value !Value of initial guess (current + continuation)

    integer :: optim_subroutine = 1 !Temporarily place here,
    !later on put into parameter file. 1 means e04jcf (works on HPC),
    !2 means e04jyf
    !
    !1 leads to crashes on HPC for some reasons (maybe older version of the NAg library)

    real(dp), dimension(ncvar_par) :: xc_perm
    real(dp), parameter :: perm_eps = 0.2_dp
    integer :: c_ind
    real(dp), dimension(M_par,I_par) :: c_guess_perm
    real(dp), dimension(ncvar_par) :: best_xc
    real(dp) :: best_fmax

    !Current_period choices
    real(dp), dimension(M_par,I_par) :: lab_mat,b_prime,a_prime,u_c,u_l
    real(dp), dimension(M_par,1) :: rho_prime,R

    real(dp), dimension(1,I_par) :: y_agent !product of agent
    real(dp), dimension(1,I_par) :: allinc_agent !all income of agent (labour + net income from assets)

    real(dp), dimension(I_par) :: a_min_slack !slackness of constraint on range of
    !a_prime from below which follows from choice of grid (pars%a_min), for all agents a_prime
    !(positive values mean that the constraint is strictly satisfied, negative mean the opposite, 0 means that the
    !constraint is satisfied with equality).
    real(dp), dimension(I_par) :: a_max_slack !the same but for constraint on a from above.
    real(dp), dimension(:), allocatable :: a_slack_all !contains the lowest slackness of all constraints
    !on agent's a_prime in a period). I.e., the most violated constraint in a period is reflected in this,
    !along with seriousness of the violation.

    integer :: att_counter

    !Common block to be accessible by subroutine CPC_ret
    common /cpc_ret_bloc/ a_st,theta_lvl,rho_st,t_st,g_st,g_st_ind
    !(EFE) : The bounds could actually be prepared even before this subroutine is called,
    !they are the same at every gridpoint.
    common /cpc_ret_bloc2/ g_lvl_mat, P_onerow,l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec

    common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
    !to some subroutines.

    optim_subroutine = pars%opt_subroutine

    !If the number of free variables in the optimization problem is 1, then
    !we can't use subroutine e04jcf. This will be the case only when we are solving
    !deterministic version of the model with two players.
    if(M_par*I_par-1 == 1) then
        optim_subroutine = 2
    end if

!_____________________________________________________________________________
    !Save the first-period results computed previously as first elements of the simulated series
    sim_result%c_sim(1,:) = FP_result%c(1,:)
    sim_result%u_c_sim(1,:) = FP_result%u_c(1,:)
    sim_result%u_l_sim(1,:) = FP_result%u_l(1,:)
    sim_result%l_sim(1,:) = FP_result%l(1,:)
    sim_result%b_prime_sim(1,:) = FP_result%b_prime(1,:)
    sim_result%a_prime_sim(1,:) = FP_result%a_prime(1,:)
    sim_result%rho_prime_sim(1) = FP_result%rho_prime
    sim_result%t_prime_sim(1) = FP_result%t_prime

    !Get t_sim and b_sim from the variables with _prime suffix at the end,
    !the difference being that these will contain the initial conditions as well.
    !Here we store some data duplicitly but it's not a problem and doing it
    !like this is the best way to avoid any confusion.
    !Theta_lvl_sim will be computed from theta_rel_sim at the end.

    !First elements of shock series - recovered from FP_result, not directly from parameter
    !file. Because if we get this from parameter file, we are going to have an
    !expected trend in g. It's better to draw this when solving first-period problem
    !from a stationary distribution.
    sim_result%s_ind_sim(1) = FP_result%s_ind_init

    sim_result%s_real_sim(1) = pars%S(1,sim_result%s_ind_sim(1))
    !Draw chain of shock realizations and save it in sim_result
    !(don't save G_lvl yet - fill these in gradually).
    !Allocate memory (we don't need T_sim but T_sim-1 elements, already have the first one)
    sim_result%s_ind_sim(2:pars%T_sim) = markov(FP_result%s_ind_init,pars%P,pars%T_sim-1)
    sim_result%s_real_sim(2:pars%T_sim) = markov_ind_to_real(reshape(pars%S,[size(pars%S)]),&
    sim_result%s_ind_sim(2:pars%T_sim))

    !Get initial guess of consumption for period t=1 (in here the corresponding index
    !will be t=2 as the elements of sim_series are indexed from 1 and not from 0).
    c_guess = matmul(ones_col,FP_result%c)
    !(If shock are large, something more sophisticated could be needed)

    !prepare bounds for checking constraints
    call prep_bounds(pars%l_max,pars%b_min,pars%b_max,pars%a_min,pars%a_max,pars%rho_min,pars%rho_max,&
    l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec)

    !In periods (indexed from t=2) to T_sim_max, already have intiial guess for c.
    !(don't load it from FP_result as only the actually realized consumption is saved
    !there, simply use c_guess from last-period.

    do t_ind=2,pars%T_sim
        c_guess_old = c_guess !back up the original guess
        !(c_guess will be the best choice found in the previous period)

        !Initialize variables related to best choice.
        best_fmax = 1000000000000000.0_dp

        att_counter = 0

167     if(att_counter == 0) then
            c_guess = c_guess_old !use the previous-period consumption plan as the initial guess
        else !permutate initial guess
            call random_number(xc_perm)
            do c_ind = 1,ncvar_par
                xc_perm(c_ind) = 1.0_dp - perm_eps + xc_perm(c_ind)*(2.0_dp*perm_eps)
            end do
            !transform into matrix form
            call c_inc_vtm(xc_perm,c_guess_perm)
            c_guess = c_guess_old * c_guess_perm
        end if

        !If we are in later attempts (this will be the case for about half of attempts and more
        !if the optimization returns subroutines with negative consumption) use fix_c (which
        !will force the initial guess to imply consumption which is independent of shock realization).
        !As in other places (such as LFFC and acc_test), this should provide a fairly good initial guess,
        !because for small shocks the variation of consumption across states should be really quite small.
        if(att_counter >= ceiling(real(pars%sim_attmax2,dp)/2.0_dp)) then
            call fix_c(c_guess,rho_st,pars)
        end if


        att_counter = att_counter + 1

        !Get states from previous elements of sim_result (so we can use the same code as elsewhere)
        a_st = sim_result%a_prime_sim(t_ind-1,:)
        rho_st = sim_result%rho_prime_sim(t_ind-1)
        t_st = sim_result%t_prime_sim(t_ind-1)
        g_st = sim_result%s_real_sim(t_ind-1) !actual shock realization (not the index)
        !(the last-period one, serving as the state this period).
        g_st_ind = sim_result%s_ind_sim(t_ind-1)

        !Get productivity in levels (the same for all what_binds and consumption choices)
        call prod_lvl(t_st,pars%theta_0,pars%xi_theta,pars%static_environment,theta_lvl)

        !Need government expenditure in levels for all shock realizations. Actually save it in a MxI matrix
        !where every row contains gov. exp. in level in state m
        do m_ind = 1,M_par
            g_lvl_mat(m_ind,:) = pars%S(1,m_ind)*trend_g(pars%mass,theta_lvl,pars%l_max,pars%k_gov)
        end do

        call max_feas_cons(pars%l_max,pars%mass,theta_lvl,g_lvl_mat,max_cons_mat)
        !Use this to generate the vector of upper bounds...

        !get upper bounds
        call c_inc_mtv(max_cons_mat,bu)

        !This should be saved in pars so can save a tiny bit of time by loading this
        !P_onerow = matmul(reshape([1.0_dp,1.0_dp],[M_par,1]),reshape(pars%P(g_st_ind,:),[1,M_par]))
        P_onerow = pars%P_onerow(g_st_ind,:,:)

        !Transform the initial guess into vector form so that it can be used
        !by NAG subroutines.
        call c_inc_mtv(c_guess,xc)

        !to make sure best_xc is properly initialized (even in a pathological case)
        if(att_counter == 1) best_xc = xc

        iuser = 0 !nothing binds, interior solution. If we want to generalize to corner solutions, then
        !we would need to cycle over what_binds as in VFI_solve_gp. (this variable
        !is used to pass what_binds to some subroutines called by optimizaiton subroutines.
        ifail = 1 !silent exit would be 1
        select case(optim_subroutine)
            case(1) !Quadratic approximation method
                call e04jcf(CPC_ret2,M_par*I_par-1,ncvar_par*2+1,xc,bl,bu,0.1_dp,&
                    pars%sim_rhoend,e04jcp,pars%sim_maxiter,fmax,nf,iuser,ruser,ifail)
            case(2) !Newton-like method
                call e04jyf(M_par*I_par-1,0,CPC_ret,bl,bu,xc,fmax,iw,liw,w,lw,iuser,ruser,ifail)
            case default
                write(*,*) 'Error in subroutine sim_series - wrong value of optim_subroutine'
                error stop
        end select
        !write(*,*) ifail
        !error stop

            !If the choice found is better than the best choice found so far, save it.
            if(fmax < best_fmax) then
                best_xc = xc
                best_fmax = fmax
            end if
            !If we haven't reached the maximum number of attempts, try again with a different initial guess
            if(att_counter < pars%sim_attmax2) go to 167

            !Take the best choice and function value
            xc = best_xc
            fmax = best_fmax

        c_im_found = .false.
        !Now we have the solution in xc. Need to get the full consumption plan and save
        !it in c_guess so it is used the next period. Assume interior solution here.
        !First get the incomplete consumption matrix (element M_par,I_par missing)
        call c_inc_vtm(xc,c_guess)
        !(iuser contains what_binds, in this case one for interior solution)
        iuser = 0 !just in case some NAg subroutine overwrites it

        call get_c_im(c_guess,iuser,&
            a_st,rho_st,theta_lvl,g_lvl_mat(:,1),g_st_ind,pars,P_onerow,c_im_found)



        !If nothing found, we have a serious problem (consumption negative for some agent)
        if(.not. c_im_found) then

            !Try up to 20 more times to get a different solution. If we still haven't found anything that does
            !not violate non-negativity, then use last-period solution.

            if(att_counter < pars%sim_attmax2 + 20) then
                go to 167
            end if

            c_guess = c_guess_old !get the thing which didn't lead to negative c_{i,M} or was the initial period guess

            write(*,*) 'Warning: get_c_im found no solutions so we cannot get a full consumption plan to be saved.'
            write(*,*) 'Check subroutine sim_series.'
            !error stop
        end if
        !Having the consumption plan (in c_guess), we need to recover plan for all variables of interest.

        call lab(c_guess,lab_mat,theta_lvl,pars%mass,g_lvl_mat(:,1),pars%Gamma_par)

        !Compute next-period states (need to compute a_prime and theta_prime)
        t_st_pr = t_prime(t_st,grids%t_gr(pars%N_t))

        call ass(a_st,c_guess,lab_mat,pars%A_par,pars%B_par,pars%Gamma_par,&
        pars%beta,P_onerow,a_prime)


        call get_b_prime(b_prime,a_prime,c_guess,pars)

        call util_c(c_guess,u_c,pars%A_par)
        rho_prime = reshape(u_c(:,1)/u_c(:,2),[M_par,1])

        !We've got plans for most of the variables which we save. The current-period
        !shock index is sim_series%s_ind_sim(t_ind). Take that row of the plans and save it in
        !c_sim_series.
        sim_result%a_prime_sim(t_ind,:) = a_prime(sim_result%s_ind_sim(t_ind),:)
        sim_result%rho_prime_sim(t_ind) = rho_prime(sim_result%s_ind_sim(t_ind),1)
        sim_result%t_prime_sim(t_ind) = t_st_pr
        !(we already have generated and saved a series of shocks)

        !We have saved the variables essential for solving the model (states)
        !Now we can compute other things and save them into sim_result
        !(government expenditure in levels, taxes, etc.)

        sim_result%c_sim(t_ind,:) = c_guess(sim_result%s_ind_sim(t_ind),:)
        sim_result%l_sim(t_ind,:) = lab_mat(sim_result%s_ind_sim(t_ind),:)
        sim_result%u_c_sim(t_ind,:) = u_c(sim_result%s_ind_sim(t_ind),:)
        sim_result%b_prime_sim(t_ind,:) = b_prime(sim_result%s_ind_sim(t_ind),:)
        !MU of labour
        call util_l(lab_mat,u_l,pars%B_par,pars%Gamma_par)
        sim_result%u_l_sim(t_ind,:) = u_l(sim_result%s_ind_sim(t_ind),:)

        !Get (gross) rate of return on the bonds
        !in periods t, when we know period t consumption plan, we compute period
        !t-1 rate of return. The length of the simulated series for R_t will this be
        !one less than the length of the other series. Conceptually we could easily
        !get the missing element but it is cumbersome.

        !cycle over shock realization to compute R (nobody cares about efficiency here)
        !can choose any agent i to compute this.
        do m_ind = 1,M_par
            R(m_ind,1) = (1.0_dp/pars%beta)*u_c(m_ind,1)/&
            maxval(matmul(pars%P(sim_result%s_ind_sim(t_ind-1):sim_result%s_ind_sim(t_ind-1),:),u_c(:,1:1)))
            !The maxval is there only to conver the (1,1) array to a scalar
        end do
        !The interest rate depends on last-period shock realization only (which
        !is the state in period t)
        sim_result%R_sim(t_ind-1) = R(sim_result%s_ind_sim(t_ind-1),1)
    end do






    !Now we have some series saved in sim_result. We can compute a few others
    sim_result%b_sim(1,:) = pars%b_init(1,:)
    sim_result%b_sim(2:pars%T_sim,:) = sim_result%b_prime_sim(1:pars%T_sim-1,:)

    sim_result%t_sim(1) = 0.0_dp !'time' variable in first period is 0
    sim_result%t_sim(2:pars%T_sim) = sim_result%t_prime_sim(1:pars%T_sim-1)

    do t_ind=1,pars%T_sim
        !Get and save productivity in levels
        call prod_lvl(sim_result%t_sim(t_ind),pars%theta_0,pars%xi_theta,&
        pars%static_environment,sim_result%theta_lvl_sim(t_ind,:))

        !sim_result%theta_lvl_sim(t_ind,:) contains productivity in level.
        sim_result%g_lvl_sim(t_ind) = trend_g(pars%mass,sim_result%theta_lvl_sim(t_ind,:)&
            ,pars%l_max,pars%k_gov)*sim_result%s_real_sim(t_ind)

        !Get the tax rate
        sim_result%tau_sim(t_ind) = sim_result%u_l_sim(t_ind,1)/(sim_result%u_c_sim(t_ind,1)&
        *sim_result%theta_lvl_sim(t_ind,1)) + 1.0_dp

        !Compute product of each agent (not saved)
        y_agent(1,:) = sim_result%theta_lvl_sim(t_ind,:)*sim_result%l_sim(t_ind,:)

        !Get Gini coefficients in terms of productivity and in terms of labour income
        !(the former is interesting only for calibration purposes, what we
        !really want to know are equilibrium outcomes).
        !sim_result%gini_income_sim,
        call gini_i2(reshape(y_agent,[I_par]),reshape(pars%mass,[I_par]),sim_result%gini_income_sim(t_ind))

        call gini_i2(reshape(sim_result%theta_lvl_sim(t_ind,:),[I_par]),reshape(pars%mass,[I_par]),&
        sim_result%gini_prod_sim(t_ind))

        !Get GDP, government assets, and Gov assets to GDP
        sim_result%Y_sim(t_ind) = sum(sim_result%theta_lvl_sim(t_ind,:) * sim_result%l_sim(t_ind,:) * &
        reshape(pars%mass,[I_par]))

        sim_result%B_gov_sim(t_ind) = sum(sim_result%b_sim(t_ind,:) * reshape(pars%mass,[I_par])) * (-1.0_dp)

        sim_result%B_gov_to_GDP_sim(t_ind) = sim_result%B_gov_sim(t_ind) / sim_result%Y_sim(t_ind)

        !Get government expenditure to GDP G_to_Y_sim
        sim_result%G_to_Y_sim(t_ind) = sim_result%g_lvl_sim(t_ind) / sim_result%Y_sim(t_ind)


        !get income of agent including assets (labour income + Rb_{t-1} -b_t)
        if(t_ind > 1) then
            allinc_agent(1,:) = sim_result%b_sim(t_ind,:)*sim_result%R_sim(t_ind-1) - sim_result%b_prime_sim(t_ind,:)
        else !R_{t-1}=1 in the first period (assumption) (could be anything really)
            allinc_agent(1,:) = sim_result%b_sim(t_ind,:) - sim_result%b_prime_sim(t_ind,:)
        end if
        !We have computed all_inc_agent for the current period only. Now use it to compute GINI
        !coefficient for all income and save it.
        !Add labour income to net asset income
        allinc_agent(1,:) = allinc_agent(1,:) + y_agent(1,:)
        call gini_i2(reshape(allinc_agent,[I_par]),reshape(pars%mass,[I_par]),sim_result%gini_allincome_sim(t_ind))

    end do


    !Now that we have the simulated series, compute a measure of 'goodness of solution. with respect
    !to hitting bounds for state a_prime (bounds for rho are not an issue and the scale is different anyhow).

    !The measure a_slack_all_avg will be average (over periods of simulation) slackness wrt constraints
    !on a_prime (a_min, a_max given in parameter file). Negative values indicate violation of constraint.
    !in every period the min slackness is computed (from slackness for every agent and both min and max constraint),
    !which determines the worst violation in period. This is saved in a_slack_all, which is then averaged.
    !The lower the average the worse the solution (hitting bounds which are arbitrary and have no economic
    !meaning).
    allocate(a_slack_all(pars%T_sim))

    !In every period, compute the slack wrt violating constraints
    do t_ind = 1,pars%T_sim
        !Compute the slackness wrt a_min for agents 1 and 2
        a_min_slack(1) = sim_result%a_prime_sim(t_ind,1) - pars%a_min(1,1)
        a_min_slack(2) = sim_result%a_prime_sim(t_ind,2) - pars%a_min(1,2)

        !Compute the slackness wrt a_max for agents 1 and 2
        a_max_slack(1) = pars%a_max(1,1) - sim_result%a_prime_sim(t_ind,1)
        a_max_slack(2) = pars%a_max(1,2) - sim_result%a_prime_sim(t_ind,2)

        !a_slack_all contains the minimum slackness from all the constraints (and thus captures
        !the biggest issue in the sense of largest violation or getting closest to the bound)
        a_slack_all(t_ind) = min(a_min_slack(1),a_min_slack(2),a_max_slack(1),a_max_slack(1))

    end do

    sim_result%a_slack_all_avg = 0.0_dp
    do t_ind = 1,pars%T_sim
        sim_result%a_slack_all_avg = sim_result%a_slack_all_avg + a_slack_all(t_ind)/real(pars%T_sim,dp)
    end do

    deallocate(a_slack_all)

end subroutine sim_series


!Subroutine sim_stat computes various statistics from multiple simulations (computed beforehand).
subroutine sim_stat(sim_res_all,pars,sim_avg)
    use nag_library, Only: G01AMF

    type(sim_result_all_type), intent(inout) :: sim_res_all !on the input contains all simulated series,
    !some additional results will be added (such as residuals, statistics, etc.)
    type(par), intent(in) :: pars
    type(sim_result_type), intent(inout) :: sim_avg !will contain the sample average of all series

    integer :: sim_ind

    real(dp) :: min_tau_res,max_tau_res !minimum and maximum residual in tau
    real(dp) :: min_tau,max_tau !min and max tax rate in the average series.

    type(par) :: pars_comm !copy of pars to be passed to subroutines for manipulating series
    !which are called using operator overloading (hence can't be passed as an input directly
    !or at least I don't know how to do that)

    real(dp), allocatable, dimension(:,:) :: tau_resid !residuals in tax policy (relative to the average realization)
        !First dimension will be the sample index, the second will be the time index in the series.

    real(dp) :: tmp_avg1,tmp_avg2 !used for computing averages which don't need to be stored...

    integer :: m_ind,n_ind,t_ind !indices for cycling over shock realizations, samples, and periods
    real(dp), allocatable, dimension(:) :: y_all_m !will contain all Y's which resulted for a given realization of shock.
    integer :: num_of_m !number of realizations of a given shock index

    real(dp), allocatable, dimension(:) :: tmp_vector !for computing quantiles

    !some variables used by NAG subroutine used to compute quantiles:
    integer :: n,nq,ifail
    real(kind = nag_wp) :: Q(2),QV(2)

    common /pars_comm/ pars_comm

    pars_comm = pars

    !First compute the sample average series.
    !(this uses operator overloading)
    sim_avg = sim_div(sim_res_all%SRA(1),real(pars%N_sim,dp))

    do sim_ind = 2,pars%N_sim
        sim_avg = sim_avg + sim_div(sim_res_all%SRA(sim_ind),real(pars%N_sim,dp))
    end do
    !The definition of operators / and * is wrong - there must be another way to define
    !binary operators for variables of different types (such as custom type and scalar
    !but it's not a priority)

    !Compute residuals (in tax policy)
    allocate(tau_resid(pars%N_sim,pars%T_sim))
    do sim_ind = 1,pars%N_sim
        tau_resid(sim_ind,:) = sim_res_all%SRA(sim_ind)%tau_sim - sim_avg%tau_sim
    end do

    !Compute the maximum absolute difference between residuals in a sample
    !For measuring within sample variation (
    do sim_ind = 1,pars%N_sim
        min_tau_res = minval(tau_resid(sim_ind,:))
        max_tau_res = maxval(tau_resid(sim_ind,:))

        sim_res_all%tau_resid_MAD(sim_ind) = abs(max_tau_res-min_tau_res)

        !Also compute variance (within sample)
        sim_res_all%tau_resid_var(sim_ind) = var_series(reshape(tau_resid(sim_ind,1:pars%T_sim),[pars%T_sim]))
    end do

    !Now compute MAD and variance of the average (across samples) series tau (measure of importance
    !of trend which is mostly driven by hetero. prod. growth). It might also be useful to take a
    !trend from deterministic environment for this (because the trend in the stochastic environment is
    !also affected by the shocks) - the only issue is that if we change the environment, the economy is a bit
    !different due to effect of levels (but hopefully it should not matter too much).
        min_tau = minval(sim_avg%tau_sim)
        max_tau = maxval(sim_avg%tau_sim)
        sim_res_all%avg_tau_MAD = abs(max_tau - min_tau)
        sim_res_all%avg_tau_var = var_series(sim_avg%tau_sim)

    !Save the residuals in the last period - this gives a sense of influence of residuals
    !because the 'residuals' in the last period capture the cumulative influence of previous
    !period shocks. Can do stuff like plot 95% confidence intervals, etc.
    do sim_ind = 1,pars%N_sim
        sim_res_all%tau_resid_last(sim_ind) = tau_resid(sim_ind,pars%T_sim)
    end do

    !The following are cross-sectional statistics on residuals in the last period.
    !If these are relatively low it means that the cumulative influence of shocks is
    !relatively low.
    !But still the influence of shocks could be reflected in the trend! So need to be careful
    !about this!!!
        min_tau_res = minval(sim_res_all%tau_resid_last)
        max_tau_res = maxval(sim_res_all%tau_resid_last)
        sim_res_all%tau_resid_CS_MAD = abs(max_tau_res - min_tau_res)
        sim_res_all%tau_resid_CS_var = var_series(sim_res_all%tau_resid_last)

    !Also save average tax in period 1, T_sim, and their difference  (tau(T_sim) - tau(1)).
    sim_res_all%tau_first = sim_avg%tau_sim(1)
    sim_res_all%tau_last = sim_avg%tau_sim(pars%T_sim)
    sim_res_all%tau_lastminfirst = sim_res_all%tau_last - sim_res_all%tau_first

    !No transfers in this model - save zeroes though, so the same script can be used for plotting as
    !for the model with transfers.
    sim_res_all%trans_first = 0.0_dp
    sim_res_all%trans_last = 0.0_dp
    sim_res_all%trans_lastminfirst = 0.0_dp


    !What follows are measures of importance of shocks
    !1) within sample (relative to trend) - in terms of MAD and var ratios
        !avg(tau_resid_MAD)/avg_tau_MAD - stochasticity metric within sample (Max Abs Diff)
        tmp_avg1 = 0.0_dp !avg_tau_MAD
        do sim_ind = 1,pars%N_sim
            tmp_avg1 = tmp_avg1 + sim_res_all%tau_resid_MAD(sim_ind)/pars%N_sim
        end do
        sim_res_all%stoch_metric_wi_MAD = tmp_avg1/sim_res_all%avg_tau_MAD

        !avg(tau_resid_var)/avg_tau_var - stochasticity metric within sample (variance)
        tmp_avg2 = 0.0_dp !avg_tau_var
        do sim_ind = 1,pars%N_sim
            tmp_avg2 = tmp_avg2 + sim_res_all%tau_resid_var(sim_ind)/pars%N_sim
        end do
        sim_res_all%stoch_metric_wi_var = tmp_avg2/sim_res_all%avg_tau_var


    !2) across samples (cross-sectional variance of tau residuals in last period relative to trend
        !variation metric
        !tau_resid_CS_MAD/avg_tau_MAD - stochasticity metric between samples (Max Abs Diff)
        sim_res_all%stoch_metric_bw_MAD = sim_res_all%tau_resid_CS_MAD/tmp_avg1
        !tau_resid_CS_var/avg_tau_var - stochasticity metric between samples (variance)
        sim_res_all%stoch_metric_bw_var = sim_res_all%tau_resid_CS_var/tmp_avg2

    deallocate(tau_resid)

    !Compute average aggregate product for each shock realization (and set it to 0.0 for shock
    !realizations which are never attained)

    allocate(y_all_m(pars%T_sim*pars%N_sim)) !max space we could need if all shocks were of the same realization
    !cycle over all shock realizations
    do m_ind = 1,M_par
        num_of_m = 0
        sim_res_all%avg_y(m_ind) = 0.0_dp !initialize

        !cycle over all samples and all periods and add the product to sim_res%y_avg(m_ind)
        do n_ind = 1,pars%N_sim
            do t_ind = 1,pars%T_sim
                !if the shock realizations for this observation has index m_ind, add it to the
                !average.
                if(sim_res_all%SRA(n_ind)%s_ind_sim(t_ind) == m_ind) then
                    sim_res_all%avg_y(m_ind) = sim_res_all%avg_y(m_ind) + sim_res_all%SRA(n_ind)%Y_sim(t_ind)
                    num_of_m = num_of_m + 1
                end if
            end do
        end do
    !divide by the number of observations to get the sample mean
    if(num_of_m > 0) then
        sim_res_all%avg_y(m_ind) = sim_res_all%avg_y(m_ind)/real(num_of_m,dp)
    else
        sim_res_all%avg_y(m_ind) = 0.0_dp
    end if
    end do
    deallocate(y_all_m)

    !Also compute quantiles for state variables - for debugging purposes (detection of cases
    !where bounds are often hit or closely approached for some shock realizations, even though
    !the average is not close to the bounds).

    !The quantiles are computed for every period (across samples)
    allocate(tmp_vector(pars%N_sim))

    !get the quantiles for a_prime
    do t_ind = 1,pars%T_sim
        !cycle over samples and save all values of a_prime into vector
        do n_ind = 1,pars%N_sim
            tmp_vector(n_ind) = sim_res_all%SRA(n_ind)%a_prime_sim(t_ind,1)
        end do

        N = pars%N_sim !number of points in the vector
        nq = 2 !number of quantiles
        q = [0.05_dp, 0.95_dp] !quantile values
        qv = 0.0_dp
        ifail = -1

        call G01AMF(N,tmp_vector,NQ,Q,QV,IFAIL)

        sim_res_all%a_prime_05(t_ind) = qv(1)
        sim_res_all%a_prime_95(t_ind) = qv(2)
    end do

    !get the quantiles for rho_prime
    do t_ind = 1,pars%T_sim
        !cycle over samples and save all values of a_prime into vector
        do n_ind = 1,pars%N_sim
            tmp_vector(n_ind) = sim_res_all%SRA(n_ind)%rho_prime_sim(t_ind)
        end do

        N = pars%N_sim !number of points in the vector
        nq = 2 !number of quantiles
        q = [0.05_dp, 0.95_dp] !quantile values
        qv = 0.0_dp
        ifail = -1

        call G01AMF(N,tmp_vector,NQ,Q,QV,IFAIL)

        sim_res_all%rho_prime_05(t_ind) = qv(1)
        sim_res_all%rho_prime_95(t_ind) = qv(2)
    end do

    deallocate(tmp_vector)

end subroutine sim_stat



subroutine save_series(sim_result,folder_name,pars)
    use ifport !so we can use system()

    type(sim_result_type) :: sim_result
    character*(*), intent(in) :: folder_name
    type(par), intent(in) :: pars

    character*256 :: shell_cmd,file_path
    integer(4) :: sys_call_result
    logical :: folder_exists

    !If folder sim in results folder does not exist, create it
    inquire(directory = 'results/'//trim(folder_name)//'/sim',exist = folder_exists)
    if (.not. folder_exists) then
        !Create folder:
        shell_cmd = 'mkdir -p results/'//trim(folder_name)//'/sim' !-p creates the directory even if directory results doesn't exist
        sys_call_result = system(shell_cmd) !this thing uses IFPORT, compiler dependent.
    end if

    !For work with Matlab, it might be better to write the elements one by one...

    !(1) Save the one-dimensional series (most of them, some are not particularly interesting)

    file_path = 'results/'//trim(folder_name)//'/sim/s_ind.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%s_ind_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/s_real.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%s_real_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/t.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%t_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/R.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%R_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/tau.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%tau_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/g_lvl.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%g_lvl_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/rho_prime.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%rho_prime_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/gini_income.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%gini_income_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/gini_prod.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%gini_prod_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/gini_allinc.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%gini_allincome_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/Y.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%Y_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/B_gov.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%B_gov_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/B_gov_to_GDP.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%B_gov_to_GDP_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/G_to_Y.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%G_to_Y_sim
    close(unit = 20)

    !(2) Save the variables whose dimension depends on the number of agents.
    file_path = 'results/'//trim(folder_name)//'/sim/c.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%c_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/l.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%l_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/b.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%b_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/a_prime.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%a_prime_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/u_c.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%u_c_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/u_l.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%u_l_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/theta.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%theta_lvl_sim
    close(unit = 20)


    !(3)Also save M_par and I_par and T_sim (for Matlab to load) in one vector : matlab_info.out
    file_path = 'results/'//trim(folder_name)//'/sim/matlab_control.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) M_par
    write(unit=20, fmt=*) I_par
    write(unit=20, fmt=*) pars%T_sim
    write(unit=20, fmt=*) pars%N_sim
    close(unit = 20)

    !If we got here without crash it means that the results were saved successfuly
    write(*,*) 'Simulated series were saved in folder results/',trim(folder_name),'.'

end subroutine save_series

!Subroutine save_stat saves some statistics about the solution (which are saved in
!type sim_series_all, along with all the simulated series).
subroutine save_stat(sim_res_all,folder_name,pars)
    type(sim_result_all_type), intent(in) :: sim_res_all
    character*(*), intent(in) :: folder_name
    type(par), intent(in) :: pars

    character*256 :: shell_cmd,file_path
    integer(4) :: sys_call_result
    logical :: folder_exists

    integer :: sim_ind

    real(dp), dimension(pars%N_sim,pars%T_sim) :: tmp_array

    !If folder sim in results folder does not exist, create it
    inquire(directory = 'results/'//trim(folder_name)//'/sim',exist = folder_exists)
    if (.not. folder_exists) then
        !Create folder:
        shell_cmd = 'mkdir -p results/'//trim(folder_name)//'/sim' !-p creates the directory even if directory results doesn't exist
        sys_call_result = system(shell_cmd) !this thing uses IFPORT, compiler dependent.
    end if

    !Save all the results in the same folder as the rest (with prefix 'stat' so it's
    !more orderly as there are going to be a lot of files in the folder).

    !The explanation of the meaning of these variables is in the definition of type sim_result_all_type

    file_path = 'results/'//trim(folder_name)//'/sim/stat_tau_resid_last.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%tau_resid_last
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_tau_resid_MAD.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%tau_resid_MAD
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_tau_resid_var.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%tau_resid_var
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_avg_tau_MAD.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%avg_tau_MAD
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_avg_tau_var.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%avg_tau_var
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_tau_resid_CS_var.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%tau_resid_CS_var
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_tau_resid_CS_MAD.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%tau_resid_CS_MAD
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_tau_first.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%tau_first
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_tau_last.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%tau_last
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_tau_lastminfirst.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%tau_lastminfirst
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_trans_first.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%trans_first
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_trans_last.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%trans_last
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_trans_lastminfirst.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%trans_lastminfirst
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_stoch_metric_wi_MAD.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%stoch_metric_wi_MAD
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_stoch_metric_wi_var.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%stoch_metric_wi_var
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_stoch_metric_bw_MAD.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%stoch_metric_bw_MAD
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_stoch_metric_bw_var.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%stoch_metric_bw_var
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_avg_y.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%avg_y
    close(unit = 20)

    !Also save Pareto weight on agent 1 ('rich') - this is useful for plotting (could be read
    !directly from parameter file but it's easier this way)
    file_path = 'results/'//trim(folder_name)//'/sim/stat_alpha1.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) pars%alpha(1,1)
    close(unit = 20)

    !Finally, save the quantiles (across simulations) for state variables (debugging purposes)
    file_path = 'results/'//trim(folder_name)//'/sim/stat_a_prime_05.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%a_prime_05
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_a_prime_95.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%a_prime_95
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_rho_prime_05.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%rho_prime_05
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_rho_prime_95.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%rho_prime_95
    close(unit = 20)


    !save the tax rate for all series
    do sim_ind = 1,pars%N_sim
        tmp_array(sim_ind,:) = sim_res_all%SRA(sim_ind)%tau_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/tau_all.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array
    close(unit = 20)

    !save the tax rate for all series
    do sim_ind = 1,pars%N_sim
        tmp_array(sim_ind,:) = sim_res_all%SRA(sim_ind)%B_gov_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/B_gov_all.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array
    close(unit = 20)

    !save the tax rate for all series
    do sim_ind = 1,pars%N_sim
        tmp_array(sim_ind,:) = sim_res_all%SRA(sim_ind)%B_gov_to_GDP_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/B_gov_to_GDP_all.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array
    close(unit = 20)

end subroutine save_stat

!Function t_prime is a trivial function capturing the law of motion of 'time' state. t is incremented by one
!unless it would exceed t_max, in which case it will be set equal to t_max.
function t_prime(t,t_max)
    real(dp) :: t_prime
    real(dp), intent(in) :: t,t_max

    t_prime = minval([t + 1.0_dp,t_max])

end function t_prime


!Subroutine det_ss computes values of some variables in deterministic steady state. This only makes
!sense in the case where static_environment == .true., or after growth has ended. However, in the
!latter case the result is conditional on the initial level of assets which we do not know.
!
!So far the subroutine assumes 0 initial asset holdings. If the asset holdings are not zero for
!at least one of the agents, the system of equations needs to be solved using an iterative method.
!
!Also it assumes utility function log in c and CRRA in leisure! Generalizations will also
!yield system of equations which can't be solved analytically.
!
!Also at this stage the algorithm assumes two players!!!
subroutine det_ss(pars,theta_lvl,b_init,g_lvl,c_ss,l_ss)
    real(dp), dimension(1,I_par), intent(in) :: theta_lvl, b_init !can be different than in parameter file
    type(par), intent(in) :: pars
    real(dp), intent(in) :: g_lvl

    real(dp), dimension(1,I_par), intent(out) :: c_ss,l_ss

    real(dp) :: K1,K2 !constants

    if(I_par/=2) then
        write(*,*) 'Error: subroutine det_ss implemented for I=2 only'
        error stop
    end if

    if(maxval(abs(b_init))>0.0_dp) then
        write(*,*) 'Error: b_init must be 0 in subroutine det_ss'
        error stop
        !generalize later if needed
    end if

    !First compute labour (for log utility it does not depend on productivity!)
    l_ss(1,1) = (pars%A_par/pars%B_par)**(1/(1+pars%gamma_par))
    l_ss(1,2) = (pars%A_par/pars%B_par)**(1/(1+pars%gamma_par))

    K1 = (theta_lvl(1,1)*((l_ss(1,2))**pars%gamma_par))/(theta_lvl(1,2)*((l_ss(1,1))**pars%gamma_par))
    K2 = pars%mass(1,1)*theta_lvl(1,1)*l_ss(1,1) + pars%mass(1,2)*theta_lvl(1,2)*l_ss(1,2) - g_lvl

    c_ss(1,2) = K2/(pars%mass(1,1)*K1 + pars%mass(1,2))
    c_ss(1,1) = K1*c_ss(1,2)

end subroutine det_ss


!Subroutine gini_i2 computes Gini coefficient for two agents.
subroutine gini_i2(x,mass,gini)
    real(dp), dimension(2), intent(in) :: x,mass
    real(dp), intent(out) :: gini

    real(dp), dimension(2) :: x_loc,mass_loc

    real(dp) :: x_mean


    x_mean = (x(1)+x(2))/2.0_dp


    if(mass(1)/=mass(2)) then
        write(*,*) 'Error: computation of gini implemented for equal mass only'
        write(*,*) 'Generalize function gini_i2'
    end if

    gini = abs(x(1)-x(2))/(2.0_dp*x_mean)
end subroutine gini_i2


subroutine loss_func(n,xc,fc,iuser,ruser)
        real (dp), intent (out) :: fc
        integer, intent (in) :: n
        real (dp), intent (inout) :: ruser(*)
        real (dp), intent (in) :: xc(n)
        integer, intent (inout) :: iuser(*) !In this instance, it is used to pass what_binds

        real(dp), dimension(:,:,:,:,:), pointer :: V_old_pntr !pointer to V_old, passed using common block

        !The following variables are passed to the function using a common block
        real(dp), dimension(:), pointer :: a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr
        real(dp), dimension(:,:), pointer :: S_pntr,P_pntr, mass_pntr
        !and transition matrix.
        real(dp), pointer :: A_par_pntr,B_par_pntr,Gamma_par_pntr
        type(par), pointer :: pars_pntr

        !Values of state variables as saved in the grid (hence the suffix st).
        real(dp), dimension(I_par) :: a_st
        real(dp), dimension(I_par-1) :: rho_st
        real(dp) :: t_st,g_st
        integer :: g_st_ind
        real(dp), dimension(1,I_par) :: theta_lvl

        real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat
        real(dp), dimension(M_par,1) :: rho_min_vec,rho_max_vec

        real(dp), dimension(M_par,I_par) :: g_lvl_mat
        real(dp), dimension(M_par,M_par) :: P_onerow

        !Local variables
        real(dp), dimension(M_par,I_par) :: c_inc_mat = 0.0_dp !incomplete consumption matrix
        !Labour
        real(dp), dimension(M_par,I_par) :: lab_mat = 0.0_dp
        !next-period states
        real(dp) :: t_st_pr
        real(dp), dimension(M_par,I_par) :: a_prime,b_prime

        logical :: c_im_found = .false., constr_ok

        integer :: g_ind !index of shock realization

        real(dp) :: loss !loss from violating constraints.

        real(dp), dimension(M_par,1) :: conval_all !continuation value in all states

        !rho_prime and u_c
        real(dp), dimension(M_par,I_par) :: u_c !marginal utility of consumption
        real(dp), dimension(M_par,1) :: rho_prime


        common /pars_pointer/ pars_pntr

        common /loss_func_bloc/ a_st,theta_lvl,rho_st,t_st,g_st,g_st_ind
        common /loss_func_bloc2/ g_lvl_mat, P_onerow,l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec

        common /comm1/ a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr,A_par_pntr,B_par_pntr,&
        Gamma_par_pntr,S_pntr,P_pntr,mass_pntr

       !Get the incomplete consumption matrix from the input vector (save it as c_inc_mat)
        call c_inc_vtm(xc,c_inc_mat)

        !Depending on what_binds (stored in iuser), compute the missing element of the consumption
        c_im_found = .false.
        call get_c_im(c_inc_mat,iuser,a_st,rho_st,theta_lvl,g_lvl_mat(:,1),g_st_ind,pars_pntr,P_onerow,c_im_found)
        !(from this point on, c_inc_mat is actually the complete consumption plan)

        if(.not. c_im_found) then
            fc = 10000000.0_dp !significant discontinuity, hopefully will almost never happen, but should guarantee
            !that such a point is never chosen.
            return
        end if

        !Compute Labour (for that we need the state passed in a common block)
        call lab(c_inc_mat,lab_mat,theta_lvl,mass_pntr,g_lvl_mat(:,1),Gamma_par_pntr)

        !Compute next-period states (need to compute a_prime and theta_prime)
        t_st_pr = t_prime(t_st,t_gr_pntr(pars_pntr%N_t))

        !(EFE) : Get rid of this and use rel_prod_prime instead, more efficient here...

        call ass(a_st,c_inc_mat,lab_mat,pars_pntr%A_par,pars_pntr%B_par,pars_pntr%Gamma_par,&
        pars_pntr%beta,P_onerow,a_prime)

        call get_b_prime(b_prime,a_prime,c_inc_mat,pars_pntr)

        constr_ok = .false.
        !Check constraints (Get a loss value). All of the matrices/vectors for checking constraints need to
        !be precomputed and passed to this thing by common block.


        !Unlike in CPC_ret and its modifications, we want to decrease loss for all constraints that
        !are satisfied. That's why the last argument of the call to check_constr is .true.

        call check_constr(c_inc_mat,lab_mat,b_prime,b_min_mat,b_max_mat,l_max_mat,theta_lvl,a_prime,a_min_mat&
        ,a_max_mat,iuser,rho_min_vec,rho_max_vec,constr_ok,loss,pars_pntr%A_par,.false.,.true.)
        !Now add discounted expected continuation value: This will be computed by function dECV

        !The 'return' here is the loss function value. We want to minimize it.
        fc = loss

end subroutine loss_func

subroutine loss_func2(n,xc,fc,iuser,ruser,inform)
        real (dp), intent (out) :: fc
        integer, intent (in) :: n
        real (dp), intent (inout) :: ruser(*)
        real (dp), intent (in) :: xc(n)
        integer, intent (inout) :: iuser(*) !In this instance, it is used to pass what_binds
        integer, intent(out) :: inform

        real(dp), dimension(:,:,:,:,:), pointer :: V_old_pntr !pointer to V_old, passed using common block

        !The following variables are passed to the function using a common block
        real(dp), dimension(:), pointer :: a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr
        real(dp), dimension(:,:), pointer :: S_pntr,P_pntr, mass_pntr
        !and transition matrix.
        real(dp), pointer :: A_par_pntr,B_par_pntr,Gamma_par_pntr
        type(par), pointer :: pars_pntr

        !Values of state variables as saved in the grid (hence the suffix st).
        real(dp), dimension(I_par) :: a_st
        real(dp), dimension(I_par-1) :: rho_st
        real(dp) :: t_st,g_st
        integer :: g_st_ind
        real(dp), dimension(1,I_par) :: theta_lvl

        real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat
        real(dp), dimension(M_par,1) :: rho_min_vec,rho_max_vec

        real(dp), dimension(M_par,I_par) :: g_lvl_mat
        real(dp), dimension(M_par,M_par) :: P_onerow

        !Local variables
        real(dp), dimension(M_par,I_par) :: c_inc_mat = 0.0_dp !incomplete consumption matrix
        !Labour
        real(dp), dimension(M_par,I_par) :: lab_mat = 0.0_dp
        !next-period states
        real(dp) :: t_st_pr
        real(dp), dimension(M_par,I_par) :: a_prime,b_prime

        logical :: c_im_found = .false., constr_ok

        integer :: g_ind !index of shock realization

        real(dp) :: loss !loss from violating constraints.

        real(dp), dimension(M_par,1) :: conval_all !continuation value in all states

        !rho_prime and u_c
        real(dp), dimension(M_par,I_par) :: u_c !marginal utility of consumption
        real(dp), dimension(M_par,1) :: rho_prime


        common /pars_pointer/ pars_pntr

        common /loss_func_bloc/ a_st,theta_lvl,rho_st,t_st,g_st,g_st_ind
        common /loss_func_bloc2/ g_lvl_mat, P_onerow,l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec

        common /comm1/ a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr,A_par_pntr,B_par_pntr,&
        Gamma_par_pntr,S_pntr,P_pntr,mass_pntr

        inform = 0

       !Get the incomplete consumption matrix from the input vector (save it as c_inc_mat)
        call c_inc_vtm(xc,c_inc_mat)

        !Depending on what_binds (stored in iuser), compute the missing element of the consumption


        c_im_found = .false.
        call get_c_im(c_inc_mat,iuser,a_st,rho_st,theta_lvl,g_lvl_mat(:,1),g_st_ind,pars_pntr,P_onerow,c_im_found)
        !(from this point on, c_inc_mat is actually the complete consumption plan)

        if(.not. c_im_found) then
            fc = 10000000.0_dp !significant discontinuity, hopefully will almost never happen, but should guarantee
            !that such a point is never chosen.
            return
        end if

        !Compute Labour (for that we need the state passed in a common block)
        call lab(c_inc_mat,lab_mat,theta_lvl,mass_pntr,g_lvl_mat(:,1),Gamma_par_pntr)

        !Compute next-period states (need to compute a_prime and theta_prime)
        t_st_pr = t_prime(t_st,t_gr_pntr(pars_pntr%N_t))

        !(EFE) : Get rid of this and use rel_prod_prime instead, more efficient here...

        call ass(a_st,c_inc_mat,lab_mat,pars_pntr%A_par,pars_pntr%B_par,pars_pntr%Gamma_par,&
        pars_pntr%beta,P_onerow,a_prime)

        call get_b_prime(b_prime,a_prime,c_inc_mat,pars_pntr)

        constr_ok = .false.
        !Check constraints (Get a loss value). All of the matrices/vectors for checking constraints need to
        !be precomputed and passed to this thing by common block.


        !Unlike in CPC_ret and its modifications, we want to decrease loss for all constraints that
        !are satisfied. That's why the last argument of the call to check_constr is .true.

        call check_constr(c_inc_mat,lab_mat,b_prime,b_min_mat,b_max_mat,l_max_mat,theta_lvl,a_prime,a_min_mat&
        ,a_max_mat,iuser,rho_min_vec,rho_max_vec,constr_ok,loss,pars_pntr%A_par,.false.,.true.)
        !Now add discounted expected continuation value: This will be computed by function dECV

        !The 'return' here is the loss function value. We want to minimize it.
        fc = loss



end subroutine loss_func2


subroutine loss_func_1d(xc,fc,iuser,ruser)
        real (dp), intent (out) :: fc
        real (dp), intent (inout) :: ruser(*)
        real (dp), intent (in) :: xc
        integer, intent (inout) :: iuser(*) !In this instance, it is used to pass what_binds

        real(dp), dimension(:,:,:,:,:), pointer :: V_old_pntr !pointer to V_old, passed using common block

        !The following variables are passed to the function using a common block
        real(dp), dimension(:), pointer :: a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr
        real(dp), dimension(:,:), pointer :: S_pntr,P_pntr, mass_pntr
        !and transition matrix.
        real(dp), pointer :: A_par_pntr,B_par_pntr,Gamma_par_pntr
        type(par), pointer :: pars_pntr

        !Values of state variables as saved in the grid (hence the suffix st).
        real(dp), dimension(I_par) :: a_st
        real(dp), dimension(I_par-1) :: rho_st
        real(dp) :: t_st,g_st
        integer :: g_st_ind
        real(dp), dimension(1,I_par) :: theta_lvl

        real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat
        real(dp), dimension(M_par,1) :: rho_min_vec,rho_max_vec

        real(dp), dimension(M_par,I_par) :: g_lvl_mat
        real(dp), dimension(M_par,M_par) :: P_onerow

        !Local variables
        real(dp), dimension(M_par,I_par) :: c_inc_mat = 0.0_dp !incomplete consumption matrix
        !Labour
        real(dp), dimension(M_par,I_par) :: lab_mat = 0.0_dp
        !next-period states
        real(dp) :: t_st_pr
        real(dp), dimension(M_par,I_par) :: a_prime,b_prime

        logical :: c_im_found = .false., constr_ok

        integer :: g_ind !index of shock realization

        real(dp) :: loss !loss from violating constraints.

        real(dp), dimension(M_par,1) :: conval_all !continuation value in all states

        !rho_prime and u_c
        real(dp), dimension(M_par,I_par) :: u_c !marginal utility of consumption
        real(dp), dimension(M_par,1) :: rho_prime


        common /pars_pointer/ pars_pntr

        common /loss_func_bloc/ a_st,theta_lvl,rho_st,t_st,g_st,g_st_ind
        common /loss_func_bloc2/ g_lvl_mat, P_onerow,l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec

        common /comm1/ a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,t_gr_pntr,A_par_pntr,B_par_pntr,&
        Gamma_par_pntr,S_pntr,P_pntr,mass_pntr

       !Get the incomplete consumption matrix from the input vector (save it as c_inc_mat)
        c_inc_mat = xc

        !Depending on what_binds (stored in iuser), compute the missing element of the consumption
        c_im_found = .false.
        call get_c_im(c_inc_mat,iuser,a_st,rho_st,theta_lvl,g_lvl_mat(:,1),g_st_ind,pars_pntr,P_onerow,c_im_found)
        !(from this point on, c_inc_mat is actually the complete consumption plan)

        if(.not. c_im_found) then
            fc = 10000000.0_dp !significant discontinuity, hopefully will almost never happen, but should guarantee
            !that such a point is never chosen.
            return
        end if

        !Compute Labour (for that we need the state passed in a common block)
        call lab(c_inc_mat,lab_mat,theta_lvl,mass_pntr,g_lvl_mat(:,1),Gamma_par_pntr)

        !Compute next-period states (need to compute a_prime and theta_prime)
        t_st_pr = t_prime(t_st,t_gr_pntr(pars_pntr%N_t))

        !(EFE) : Get rid of this and use rel_prod_prime instead, more efficient here...

        call ass(a_st,c_inc_mat,lab_mat,pars_pntr%A_par,pars_pntr%B_par,pars_pntr%Gamma_par,&
        pars_pntr%beta,P_onerow,a_prime)

        call get_b_prime(b_prime,a_prime,c_inc_mat,pars_pntr)

        constr_ok = .false.
        !Check constraints (Get a loss value). All of the matrices/vectors for checking constraints need to
        !be precomputed and passed to this thing by common block.


        !Unlike in CPC_ret and its modifications, we want to decrease loss for all constraints that
        !are satisfied. That's why the last argument of the call to check_constr is .true.

        call check_constr(c_inc_mat,lab_mat,b_prime,b_min_mat,b_max_mat,l_max_mat,theta_lvl,a_prime,a_min_mat&
        ,a_max_mat,iuser,rho_min_vec,rho_max_vec,constr_ok,loss,pars_pntr%A_par,.false.,.true.)
        !Now add discounted expected continuation value: This will be computed by function dECV

        !The 'return' here is the loss function value. We want to minimize it.
        fc = loss

end subroutine loss_func_1d

!The subroutine acc_test performs accuracy tests. It is only called on one image (1).
subroutine acc_test(pars,grids,folder_name,N_t,sim_res_all)
    Use nag_library, Only: e04jyf, nag_wp, e04jcf, e04jcp,e05jaf,e05jbf,E05JBK,G01AMF!,C05AUF

    type(par), intent(in) :: pars
    type(grids_type),intent(in) :: grids
    character*(*) :: folder_name
    integer, intent(in) :: N_t

    type(sim_result_all_type), intent(in) :: sim_res_all

    real(dp), dimension(M_par,I_par) :: c_guess,c,u_c,c_bac,c_guess_perm = 1.0_dp
    real(dp), dimension(1,2) :: acc_a_min, acc_a_max
    real(dp) :: acc_rho_min, acc_rho_max
    real(dp), dimension(M_par,1) :: rho_prime, R

    integer :: samp_ind !for cycling over sample

    !Values of state variables as saved in the grid (hence the suffix st).
    real(dp), dimension(I_par) :: a_st
    real(dp), dimension(I_par-1) :: rho_st
    real(dp) :: t_st,g_st
    integer :: g_st_ind

    !prepared bounds...
    real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat
    real(dp), dimension(M_par,1) :: rho_min_vec,rho_max_vec

    real(dp), dimension(1,I_par) :: theta_lvl
    real(dp), dimension(M_par,I_par) ::  g_lvl_mat,max_cons_mat
    real(dp) :: t_st_pr
    real(dp), dimension(M_par,M_par) :: P_onerow

    integer :: optim_subroutine = 1
    integer :: m_ind

    real(dp), dimension(M_par,I_par) :: lab_mat = 0.0_dp
    !next-period states
    real(dp), dimension(M_par,I_par) :: a_prime,b_prime

    integer :: att_max,att_counter,att_max2

    character(256) :: file_path

    !Variables used in NAG subroutine E04JYF
    real(dp), dimension(ncvar_par) :: xc !choice variables (M_par x I_par matrix in row form)
    real(dp) :: fc !function value in maximization
    !The following two arrays are arrays of reals and integers which can be used to supply information
    !to the maximized function as an alternative to COMMON block.
    real(dp), dimension(1) :: ruser !Contains no information, only defined for purposes of using a NAG subroutine
    integer, dimension(3) :: iuser !will be used to pass what_binds to some subroutines.
    integer :: ifail
    !Bounds:
    real(dp), dimension(M_par*I_par-1) :: bl = 0.0001_dp,bu
    real(dp) :: fmax
    !some other things needed by  e04jyf
    !Try changing the dimension of these things and see what happens (there is only a lower bound)
    integer :: LIW = M_par*I_par +1,LW=max(ncvar_par*(ncvar_par-1)/2+12*ncvar_par,13)
    integer, dimension(M_par*I_par +1) :: IW
    real(dp), dimension(max(ncvar_par*(ncvar_par-1)/2+12*ncvar_par,13)) :: W
    integer :: nf
    logical :: c_im_found

    integer :: s_ind!, m_ind

    real(dp), dimension(ncvar_par) :: best_xc
    real(dp) :: best_fmax

    real(dp), dimension(ncvar_par) :: xc_perm
    real(dp), parameter :: perm_eps = 0.2_dp
    integer :: c_ind

    real(dp) :: rnd_num

    !period 1 and period 2 copies of the allocations.
    real(dp), dimension(M_par,I_par) :: c_t1, lab_mat_t1, a_prime_t1,g_lvl_mat_t1
    real(dp), dimension(I_par) :: a_st_t1
    real(dp), dimension(I_par-1) :: rho_st_t1
    real(dp) :: t_st_t1,g_st_t1,t_st_pr_t1
    integer :: g_st_ind_t1
    real(dp), dimension(1,I_par) :: theta_lvl_t1
    real(dp), dimension(M_par,1) :: rho_prime_t1, R_t1

    real(dp), dimension(M_par) :: omega_t1,xi_t1
    real(dp), dimension(M_par,I_par) :: lambda_e_t1


    real(dp), dimension(M_par,I_par) :: c_t2, lab_mat_t2, a_prime_t2,g_lvl_mat_t2
    real(dp), dimension(I_par) :: a_st_t2
    real(dp), dimension(I_par-1) :: rho_st_t2
    real(dp) :: t_st_t2,g_st_t2,t_st_pr_t2
    integer :: g_st_ind_t2
    real(dp), dimension(1,I_par) :: theta_lvl_t2
    real(dp), dimension(M_par,1) :: rho_prime_t2, R_t2

    !These will never be used but we need to use tham as na argument of get_lambda anyway
    real(dp), dimension(M_par) :: omega_t2,xi_t2
    real(dp), dimension(M_par,I_par) :: lambda_e_t2


    !Variables with suffix d ('delta') correspond to a solution in which consumption
    !was changed (to establish to what change in consumption the error roughly corresponds).
    real(dp), dimension(M_par,I_par) :: c_d, lab_mat_d, a_prime_d,g_lvl_mat_d
    real(dp), dimension(I_par) :: a_st_d
    real(dp), dimension(I_par-1) :: rho_st_d
    real(dp) :: t_st_d,g_st_d,t_st_pr_d
    integer :: g_st_ind_d
    real(dp), dimension(1,I_par) :: theta_lvl_d
    real(dp), dimension(M_par,1) :: rho_prime_d, R_d
    real(dp), dimension(M_par,I_par) :: lambda_d

    !multipliers in period 1 (for the actually realized shock realizations)
    real(dp), dimension(M_par,I_par) :: lambda_t1
    !and in period 2 (the last index again denotes realization of shock in period 1)
    real(dp), dimension(M_par,I_par,M_par) :: lambda_t2
    !implied lambda1
    real(dp), dimension(M_par,I_par) :: lambda_t1_imp

    real(dp), dimension(M_par,I_par) :: resid_pl !'plan for residuals' - for all states and agents
    real(dp), allocatable, dimension(:) :: resid_all !all residuals (one per initial point draw in state space,
    !at every point we average over states and agents weighted by population frequency)
    real(dp) :: resid_avg !sample average RR residual
    !The same but for individual agents (summation happens across states only). This makes sense because
    !if the level of consumption is vastly different (as is often the case) then the scale of
    !lagrange multipliers (and the errors) may be quite different...
    real(dp), allocatable, dimension(:,:) :: resid_all_ag
    real(dp), allocatable, dimension(:,:) :: resid_c_all_ag !same as above but in terms of consumption.
    real(dp), dimension(I_par) :: resid_avg_ag
    real(dp) :: resid_c1_avg !residual in terms of consumption of agent 1 (there are really upper bounds).

    !some variables used by NAG subroutine used to compute quantiles:
    integer :: n,nq
    real(kind = nag_wp) :: Q(12),QV(12)

    !Pointer to value function. This will be used to make the old value function V_old
    !accessible to the subroutine used in NAG maximization.
    real(dp), dimension(:,:,:,:,:), pointer :: V_old_pntr

    integer :: ok_counter !counts the points in state space for which no errors occured (mainly those
    !related to negative consumption) - we discard such points at the outset. These are called by bad
    !initial guess and not bad value function...

    logical :: sing_fail

    integer :: i,m
    real(dp) :: lambda_e !expectation of lambda at time t+1, contingent on s_t realization for some agent

    character(256) :: string1, string2, string3

    !Some variables needed by subroutine C05AUF
    real(dp) :: X_c, H_c, EPS_c, ETA_c, A_c, B_c

    !The following variables are passed to function c1_imp_err using a common block. They have suffix 'a' so we don't
    !mix these up with other variables
    real(dp), dimension(M_par,I_par) :: c_pl_a !consumption plan (for all agents at all states)
    integer :: s_ind_a !index of current-period shock realization.
    integer :: s_st_ind_a !index of last-period shock realization.
    real(dp) :: A_a,B_a,C_a !three pre-computed terms
    real(dp), dimension(M_par,M_par) :: P_a !probability transition matrix

    !testing RNG stuff
    integer :: rng_seed_size = 10
    integer, dimension(10) :: rng_seed
    !end of testing variables

    integer :: sim_n_ind,sim_t_ind
    real(dp) :: ind_real

    real(dp) :: best_t1
    real(dp), dimension(M_par) :: best_t2

    common /c1_imp_err_comm/ c_pl_a,A_a,B_a,C_a,P_a,s_st_ind_a,s_ind_a


    !Common block to be accessible by subroutine CPC_ret
    common /cpc_ret_bloc/ a_st,theta_lvl,rho_st,t_st,g_st,g_st_ind
    !(EFE) : The bounds could actually be prepared even before this subroutine is called,
    !they are the same at every gridpoint.
    common /cpc_ret_bloc2/ g_lvl_mat, P_onerow,l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec
    common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
    !to some subroutines.


    optim_subroutine = pars%opt_subroutine

    if(M_par*I_par-1 == 1) then
        optim_subroutine = 2
    end if

    !set number of attempts depending on which optimization routine we use
    att_max = pars%acc_attmax1
    att_max2 = pars%acc_attmax2 !in the second period the initial guess should be pretty good already so we should not need so many permutations
    !but experiment with implications of att_max and att_max2 for accuracy tests results.

    select case(pars%acc_method)
        case(1) !Euler equation residuals
            call double_output(folder_name,'Performing accuracy tests (EE residuals).')
        case default
            write(*,*) 'Error in subroutine acc_test - wrong value of acc_method in parameter file.'
            error stop
    end select

    !First we need to generate intervals for drawing samples in state space
    acc_a_min =  (pars%a_min + pars%a_max)/2.0_dp - pars%acc_int * (pars%a_max - pars%a_min)
    acc_a_max = (pars%a_min + pars%a_max)/2.0_dp + pars%acc_int * (pars%a_max - pars%a_min)
    acc_rho_min = (pars%rho_min + pars%rho_max)/2.0_dp - pars%acc_int * (pars%rho_max - pars%rho_min)
    acc_rho_max = (pars%rho_min + pars%rho_max)/2.0_dp + pars%acc_int * (pars%rho_max - pars%rho_min)

    ok_counter = 0 !initialize counter of points in state space for which we are able to obtain
    !the residuals without encountering serious errors unrelated to quality of the solution (value function).
    !The main source of these is negativity of consumption...

    allocate(resid_all(pars%acc_samples))
    allocate(resid_all_ag(I_par,pars%acc_samples))
    allocate(resid_c_all_ag(I_par,pars%acc_samples))
    resid_c_all_ag = 0.0_dp
    !allocate(EE_c1_all(pars%acc_samples))

    do samp_ind = 1,pars%acc_samples
        att_counter = 0 !initialize counter of attempts
        !initialize best_fmax at a huge value. We want to minimize this, so every time we find a choice
        !better than xc, save that choice.
        best_fmax = 1000000000000000.0_dp


        !Randomly choose a simulated series (n_ind) and then an element (a time) in the series. Use this
        !point as the state where we compute the residuals. We could also go over all of the
        !series and their elements (actually, this might be better).
        !call these indices sim_n_ind,sim_t_ind

        call random_number(rnd_num) !rnd_num is from U(0,1)
        ind_real = 1.0_dp + rnd_num*(real(pars%N_sim,dp) - 1.0_dp)
        !round to the closest real
        if(ind_real - real(floor(ind_real),dp) > 0.5_dp) then
            sim_n_ind = real(ceiling(ind_real),dp)
        else
            sim_n_ind = real(floor(ind_real),dp)
        end if


        call random_number(rnd_num) !rnd_num is from U(0,1)
        ind_real = 1.0_dp + rnd_num*(real(pars%T_sim,dp) - 2.0_dp)
        !-2.0 here, not 1.0, because we do want to assess accuracy during growth period,
        !not after growth has stopped.
        !round to the closest real
        if(ind_real - real(floor(ind_real),dp) > 0.5_dp) then
            sim_t_ind = real(ceiling(ind_real),dp)
        else
            sim_t_ind = real(floor(ind_real),dp)
        end if

        write(*,*) 'Test - imposing sim_t_ind = 1 in subroutine acc_test'
        sim_t_ind = 1

        !Get states from previous elements of sim_result (so we can use the same code as elsewhere)
        a_st = sim_res_all%SRA(sim_n_ind)%a_prime_sim(sim_t_ind,:)
        rho_st = sim_res_all%SRA(sim_n_ind)%rho_prime_sim(sim_t_ind)
        t_st = sim_res_all%SRA(sim_n_ind)%t_prime_sim(sim_t_ind)
        g_st = sim_res_all%SRA(sim_n_ind)%s_real_sim(sim_t_ind) !actual shock realization (not the index)
        !(the last-period one, serving as the state this period).
        g_st_ind = sim_res_all%SRA(sim_n_ind)%s_ind_sim(sim_t_ind)


        !Now we have the point in state space. We need to compute solution for the next period and use it
        !to compute residuals in all equations of interest
        !First do some preparatory computations
        !Get productivity in levels (the same for all what_binds and consumption choices)
        call prod_lvl(t_st,pars%theta_0,pars%xi_theta,pars%static_environment,theta_lvl)
        !Need government expenditure in levels for all shock realizations. Actually save it in a MxI matrix
        !where every row contains gov. exp. in level in state m

        !write(*,*) 'test - temporarily imposing rho_st = theta2/theta1 in acc_test'
        !rho_st = theta_lvl(1,2)/theta_lvl(1,1)
        !write(*,*) rho_st

        do m_ind = 1,M_par
            g_lvl_mat(m_ind,:) = pars%S(1,m_ind)*trend_g(pars%mass,theta_lvl,pars%l_max,pars%k_gov)
        end do
        call max_feas_cons(pars%l_max,pars%mass,theta_lvl,g_lvl_mat,max_cons_mat)
        !Use this to generate the vector of upper bounds...
        !get upper bounds
        call c_inc_mtv(max_cons_mat,bu)
        !This should be saved in pars so can save a tiny bit of time by loading this
        !P_onerow = matmul(reshape([1.0_dp,1.0_dp],[M_par,1]),reshape(pars%P(g_st_ind,:),[1,M_par]))
        P_onerow = pars%P_onerow(g_st_ind,:,:)

        !prepare bounds for checking constraints
        call prep_bounds(pars%l_max,pars%b_min,pars%b_max,pars%a_min,pars%a_max,pars%rho_min,pars%rho_max,&
        l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec)

        !Generate random initial choice of consumption.
        166 call random_number(c_guess)

        att_counter = att_counter + 1 !increase the counter of attempts made

        c_guess = c_guess*max_cons_mat*0.4_dp !the 0.4 scaling is used because for values close to max_cons_mat, those
        !choice will not be feasible anyway so it's a waste of time to try these
        !c_gp contains the actual consumption choices by agents (corresponding to c's in the paper).

        !Check whether the consumption guess implies negative consumption for agent I, state M. If yes, then
        !we do not want to start there - there is no quarantee that the algorithm will converge anywhere
        !reasonable and it's most likely going to be a huge waste of time (maximum number of iterations reached
        !in the optimization subroutines).
        call get_c_im(c_guess,iuser,a_st,rho_st,theta_lvl,g_lvl_mat(:,1),g_st_ind,pars,P_onerow,c_im_found)
        if(c_im_found == .false. .or. .true.) then
        !In the first period it is actually usually better to use fix_c (it increased accuracy a lot)
        !But only do this automatically in the first roughly half of attempts (because if the
        !number of attempts increases sufficiently it makes little sense to restrict ourselves to the
        !fix_c guesses always, as we neglect choices which are possibly better)
            call fix_c(c_guess,rho_st,pars)
        end if
        !It might actually be better to do this with all guesses - so in effect randomization would be
        !over c1 only and c2 would be chosen constant across states. This might be a much better guess
        !than a completely random one! Same in LFFC. Check the implication of doing this...

        !Transform the initial guess into vector form so that it can be used
        !by NAG subroutines.
        call c_inc_mtv(c_guess,xc)
        if(att_counter == 1) best_xc = xc
        !Just to make sure we always have something saved...

        iuser = 0 !nothing binds, interior solution. If we want to generalize to corner solutions, then
        !we would need to cycle over what_binds as in VFI_solve_gp. (this variable
        !is used to pass what_binds to some subroutines called by optimizaiton subroutines.

        ifail = 1 !silent exit would be 1 (errors do not stop execution and are not displayed).
        !!debug
!        write(*,*) 'Attempt number ',att_counter
!        ifail = -1 !errors don't stop execution but are displayed.
        !end debug

        select case(optim_subroutine)
            case(1) !Quadratic approximation method
                call e04jcf(CPC_ret2,M_par*I_par-1,ncvar_par*2+1,xc,bl,bu,0.1_dp,&
                    pars%rhoend,e04jcp,pars%acc_maxiter,fmax,nf,iuser,ruser,ifail)
            case(2) !Newton-like method
                call e04jyf(M_par*I_par-1,0,CPC_ret,bl,bu,xc,fmax,iw,liw,w,lw,iuser,ruser,ifail)
            case default
                write(*,*) 'Error in subroutine sim_series - wrong value of optim_subroutine'
                error stop
        end select


        !If the choice found is better than the best choice found so far, save it.
        if(fmax < best_fmax) then
            best_xc = xc
            best_fmax = fmax
        end if
        !If we haven't reached the maximum number of attempts, try again with a different initial guess
        if(att_counter < att_max) go to 166

        !Take the best choice and function value
        xc = best_xc
        fmax = best_fmax

        !For debugging purposes, save the best value in period 1
        best_t1 = best_fmax


        !Now we have the solution in terms of c(xc). Need to get all the other variables (just like in sim-series
        !of vfi_solve_gp) to plug into EE residual equation.

        !c is the consumption plan - incomplete at this stage
        call c_inc_vtm(xc,c)


        call get_c_im(c,iuser,a_st,rho_st,theta_lvl,g_lvl_mat(:,1),g_st_ind,pars,P_onerow,c_im_found)
        !(from this point on, c_inc_mat is actually the complete consumption plan)
        if(.not. c_im_found) then
            !skip to the next point -
            !call fix_c(c,rho_st,pars) !Consumption was probably slightly negative here.
            write(string1,*) 'Negative consumption in accuracy test in period 1..'
            call double_output(folder_name,trim(string1))
            go to 333
            !If this happens then we were really unlucky with the guesses (so the solution would be to
            !start again with random initial guesses).
            !Actually if we are not able to find a guess with non-negative consumption - just do not use this point
            !to compute residuals - it's a failure of being able to find a good initial guess, not a
            !sign that the value function is bad.
        end if
        !Compute Labour (for that we need the state passed in a common block)
        call lab(c,lab_mat,theta_lvl,pars%mass,g_lvl_mat(:,1),pars%gamma_par)
        !Compute next-period states (need to compute a_prime and theta_prime)
        t_st_pr = t_prime(t_st,grids%t_gr(pars%N_t))
        call ass(a_st,c,lab_mat,pars%A_par,pars%B_par,pars%Gamma_par,&
        pars%beta,P_onerow,a_prime)
        call get_b_prime(b_prime,a_prime,c,pars) !don't actually need b_prime but could be interesting

        !Also get other variables (rho we need as state), R we don't really need directly but it will
        !simplify things later.
        call util_c(c,u_c,pars%A_par)
        rho_prime = reshape(u_c(:,1)/u_c(:,2),[M_par,1])
        do m_ind = 1,M_par
            R(m_ind,1) = (1.0_dp/pars%beta)*u_c(m_ind,1)/&
            maxval(matmul(pars%P(g_st_ind:g_st_ind,:),u_c(:,1:1)))
            !The maxval is there only to convert the (1,1) array to a scalar
        end do

        !Save the solution into vector with suffix _t1 (time t= 1). This is a bit messy but we are using common block to
        !pass some information to CPC_ret and related subroutines so we might want to overwrite some of the
        !things. The 't+1' period solutions will be saved with suffix _t2, and they will have one more index
        !(the last index will contain the index of realization of shock in period 1, i.e., it will tell us
        !to which branch of the tree the t+1 state-contingent plan corresponds). The second reason
        !to save a copy is that we don't need to change many of the formulas below.

        c_t1 = c
        lab_mat_t1 = lab_mat
        a_prime_t1 = a_prime
        t_st_pr_t1 = t_st_pr
        theta_lvl_t1 = theta_lvl
        g_lvl_mat_t1 = g_lvl_mat
        rho_prime_t1 = rho_prime
        R_t1 = R
        a_st_t1 = a_st
        rho_st_t1 = rho_st
        t_st_t1 = t_st
        g_st_t1 = g_st
        g_st_ind_t1 = g_st_ind

        call get_lambda(lambda_t1,pars,a_st_t1,theta_lvl_t1,g_st_ind_t1,c_t1,lab_mat_t1,a_prime_t1,R_t1,sing_fail,omega_t1,xi_t1,lambda_e_t1)
        if(sing_fail) then
            write(string1,*) 'Some singular value in get_lambda less than 0.00001. Point not used.'
            call double_output(folder_name,trim(string1))
            go to 333
        end if


        !Now we need to compute allocation for all possible histories in the next period. For M=2,
        !these are two state-contingent plans.
        !Cycle over all possible shock realizations in period t. For every realization take the allocation
        !in period t as state, and compute the contingent plan for t+1.
        do g_st_ind = 1,M_par !g_st_ind is now the index of shock in period 1 (which is the state in period 2)
            !In the following periods we use the previous-period consumption plan as the initial guess.
            !Hopefully it's a good one unlike in the first period. But we still try a few permutations of
            !the initial guess and pick the best result.

            !Get state variables for this shock realization
            a_st = a_prime_t1(g_st_ind,:)
            rho_st = rho_prime_t1(g_st_ind,1)
            t_st = t_st_pr_t1 !'time' index
            g_st = pars%S(1,g_st_ind)

            call prod_lvl(t_st,pars%theta_0,pars%xi_theta,pars%static_environment,theta_lvl)

            do m_ind = 1,M_par
                g_lvl_mat(m_ind,:) = pars%S(1,m_ind)*trend_g(pars%mass,theta_lvl,pars%l_max,pars%k_gov)
            end do
            call max_feas_cons(pars%l_max,pars%mass,theta_lvl,g_lvl_mat,max_cons_mat)
            !Use this to generate the vector of upper bounds...
            !get upper bounds
            call c_inc_mtv(max_cons_mat,bu)
            P_onerow = pars%P_onerow(g_st_ind,:,:)

            !prepare bounds for checking constraints
            call prep_bounds(pars%l_max,pars%b_min,pars%b_max,pars%a_min,pars%a_max,pars%rho_min,pars%rho_max,&
            l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec)

            att_counter = 0 !initialize counter of attempts
            !initialize best_fmax at a huge value. We want to minimize this, so every time we find a choice
            !better than xc, save that choice.

            !(IMP) We could possibly compute the value of initial guess at first but if we diverge from it and pick a
            !bad value every single time, then we've got issues anyway and picking the initial guess wouldn't help much.

            best_fmax = 1000000000000000.0_dp

            167 if(att_counter == 0) then
                c_guess = c_t1 !use the first-period consumption plan as the initial guess for the 2nd period cons. plan
            else !permutate initial guess
                call random_number(xc_perm)
                do c_ind = 1,ncvar_par
                    xc_perm(c_ind) = 1.0_dp - perm_eps + xc_perm(c_ind)*(2.0_dp*perm_eps)
                end do
                !transform into matrix form
                call c_inc_vtm(xc_perm,c_guess_perm)
                c_guess = c_t1 * c_guess_perm
            end if

            att_counter = att_counter + 1

            call get_c_im(c_guess,iuser,a_st,rho_st,theta_lvl,g_lvl_mat(:,1),g_st_ind,pars,P_onerow,c_im_found)
            if(c_im_found == .false.) then
                call fix_c(c_guess,rho_st,pars)
                !In the second period using fix_c always did not improve results (and usually made them
                !slightly worse).
            end if

            !Transform the initial guess into vector form so that it can be used by NAG subroutines.
            call c_inc_mtv(c_guess,xc)

            !to make sure best_xc is properly initialized (even in a pathological case)
            if(att_counter == 1) best_xc = xc

            iuser = 0 !nothing binds, interior solution. If we want to generalize to corner solutions, then
            !we would need to cycle over what_binds as in VFI_solve_gp. (this variable
            !is used to pass what_binds to some subroutines called by optimizaiton subroutines.
            ifail = 1 !silent exit would be 1
            select case(optim_subroutine)
                case(1) !Quadratic approximation method
                    call e04jcf(CPC_ret2,M_par*I_par-1,ncvar_par*2+1,xc,bl,bu,0.1_dp,&
                        0.000001_dp,e04jcp,500,fmax,nf,iuser,ruser,ifail)
                case(2) !Newton-like method
                    call e04jyf(M_par*I_par-1,0,CPC_ret,bl,bu,xc,fmax,iw,liw,w,lw,iuser,ruser,ifail)
                case default
                    write(*,*) 'Error in subroutine sim_series - wrong value of optim_subroutine'
                    error stop
            end select

            !If the choice found is better than the best choice found so far, save it.
            if(fmax < best_fmax) then
                best_xc = xc
                best_fmax = fmax
            end if
            !If we haven't reached the maximum number of attempts, try again with a different initial guess
            if(att_counter < att_max2) go to 167

            !Take the best choice and function value
            xc = best_xc
            fmax = best_fmax

            !For debugging purposes, save the best value in period 1
            best_t2(g_st_ind) = best_fmax

            !c is the consumption plan - incomplete at this stage
            call c_inc_vtm(xc,c)
            call get_c_im(c,iuser,a_st,rho_st,theta_lvl,g_lvl_mat(:,1),g_st_ind,pars,P_onerow,c_im_found)
            !(from this point on, c_inc_mat is actually the complete consumption plan)
            if(.not. c_im_found) then
                c(M_par,I_par) = 0.0000001_dp
                write(string1,*) 'negative consumption in accuracy test in period 2.'
                call double_output(folder_name,trim(string1))
                go to 333
            end if
            !Compute Labour (for that we need the state passed in a common block)
            call lab(c,lab_mat,theta_lvl,pars%mass,g_lvl_mat(:,1),pars%gamma_par)
            !Compute next-period states (need to compute a_prime and theta_prime)
            t_st_pr = t_prime(t_st,grids%t_gr(pars%N_t))
            call ass(a_st,c,lab_mat,pars%A_par,pars%B_par,pars%Gamma_par,&
            pars%beta,P_onerow,a_prime)
            call get_b_prime(b_prime,a_prime,c,pars) !don't actually need b_prime but could be interesting

            !Also get other variables (rho we need as state), R we don't really need directly but it will
            !simplify things later.
            call util_c(c,u_c,pars%A_par)
            rho_prime = reshape(u_c(:,1)/u_c(:,2),[M_par,1])
            do m_ind = 1,M_par
                R(m_ind,1) = (1.0_dp/pars%beta)*u_c(m_ind,1)/&
                maxval(matmul(pars%P(g_st_ind:g_st_ind,:),u_c(:,1:1)))
                !The maxval is there only to convert the (1,1) array to a scalar
            end do

            !Save the results (just like in period 2 but now with suffix t2 and one more index)

            c_t2 = c
            lab_mat_t2 = lab_mat
            a_prime_t2 = a_prime
            t_st_pr_t2 = t_st_pr
            theta_lvl_t2 = theta_lvl
            g_lvl_mat_t2 = g_lvl_mat
            rho_prime_t2 = rho_prime
            R_t2 = R
            a_st_t2 = a_st
            rho_st_t2 = rho_st
            t_st_t2 = t_st
            g_st_t2 = g_st
            g_st_ind_t2 = g_st_ind

            !Actually we only need to save the lambdas for all shock realizations, not all allocation...
            call get_lambda(lambda_t2(:,:,g_st_ind),pars,a_st_t2,theta_lvl_t2,g_st_ind_t2,c_t2,lab_mat_t2,a_prime_t2,R_t2,sing_fail,omega_t2,xi_t2,lambda_e_t2)
            if(sing_fail) then
                write(string1,*) 'Some singular value in get_lambda less than 0.00001. Point not used.'
                call double_output(folder_name,trim(string1))
                go to 333
            end if

        end do

        !If we got here we successfuly computed the residual.
        !Technically not yet but in the remainder we're just computing sums of terms and there
        !is not option to fail, really.
        ! That mainly means that for the state
        !space point we got a solution with non-negative consumption in state M for agent I, and no
        !further error happened. Increase the counter of the number of solutions actually used, and save the
        !solutions into the corresponding vector (so that we may later compute averages or other
        !sample statistics).
        ok_counter = ok_counter + 1

        resid_all(ok_counter) = 0.0_dp
        resid_all_ag(:,ok_counter) = 0.0_dp
        !compute residual for all agents and all shock realizations (in period 't')
        do i = 1,I_par
        do m = 1,M_par
            !Compute expected lambda_(t+1) for agent i. This is contingent on time 't' shock realization.
            lambda_e = 0.0_dp
            do s_ind = 1,M_par
                !m is the index of period t shock realization. The t+1 transition probabilities
                !are contingent on this.
                lambda_e = lambda_e + lambda_t2(s_ind,i,m)*pars%P(m,s_ind)
            end do
            !Also get the implied lambda_t1
            lambda_t1_imp(m,i) = pars%beta * R_t1(m,1) * lambda_e

            !residual for agent i if shock realization is m
            !resid_pl(m,i) = lambda_t1(m,i) - pars%beta * R_t1(m,1) * lambda_e
            resid_pl(m,i) = lambda_t1(m,i) - lambda_t1_imp(m,i)

            !add the component of the residual 'plan' to resid_all - the saved bit is an expected
            !residual (across agents and shock realizations in period t) -conditional on g_st_ind_t1 (initially
            !drawn shock index in the state space).
            resid_all(ok_counter) = resid_all(ok_counter) + pars%mass(1,i) * pars%P(g_st_ind_t1,m) * abs(resid_pl(m,i))

            !And add to the residuals for agents only (this was added later and it could be used to compute resid_all
            !so we do not compute the expectation twice but this is pretty much costless compared to other things
            !we do in this subroutine).
            resid_all_ag(i,ok_counter) = resid_all_ag(i,ok_counter) +  pars%P(g_st_ind_t1,m) * abs(resid_pl(m,i))
        end do
        end do


        !Ok, we've got EE residuals. But we would like to get an interpretation in terms of consumption.
        !First step is to compute implied lambda_t1 - what lambda_t1 should be from Euler equations, if period
        !t+1 allocations are optimal. Then we compute what implied consumption corresponds to this.
        !(m is t-period shock realization)
        do m = 1,M_par
            !the first argument is what we vary when we are trying to find a zero. Debug - when this is c_t1(m,1),
            !and when lambda_t1_imp = lambda_t1, we should get a zero! If not, we made an error in implementing the formula.

            !We need to correctly assign/precompute the following variables to be passed using the common block.
            !c_pl_a,A_a,B_a,C_a,P_a,s_st_ind_a,s_ind_a
            c_pl_a = c_t1
            s_st_ind_a = g_st_ind_t1
            s_ind_a = m
            P_a = pars%P

            A_a = pars%mass(1,1)*pars%alpha(1,1) !ok
            B_a = lambda_t1_imp(m,1)*a_st_t1(1)*pars%P(g_st_ind_t1,m) !ok
            C_a = lambda_t1_imp(m,1)*(pars%B_par * lab_mat_t1(m,1)**(1.0_dp + pars%gamma_par) - 1.0_dp) !ok
            !The following line is different in sign for agent 2 (if we ever want to do the test for changes in agent 2's consumption)
            C_a = C_a - omega_t1(m)*( (1.0_dp/theta_lvl_t1(1,1)) * pars%B_par * lab_mat_t1(m,1)**(pars%gamma_par) ) !ok
            C_a = C_a - xi_t1(m)*pars%mass(1,1) !ok
            C_a = C_a - pars%beta**2.0_dp * R_t1(m,1) * a_prime_t1(m,1) * lambda_e_t1(m,1) !ok

            !Note: So far we're only doing the consumption interpretation for residuals in agent 1's equation.

            !We are only looking for consumption choices in the non-negative interva;
            ifail = 1 !silent exit would be 1
            !X_c, H_c, EPS_c, ETA_c, A_c, B_c
            X_c = c_t1(m,1) !initial guess - the actual consumption in period t (hopefully the implied is not too far)
            !Step length - the function will explore interval (x-256H_C,x+256H_c). We don't want to try negative consumption,
            !so we are going to try H_c = c_t1(m/1)/257
            H_c = (c_t1(m,1)*0.5_dp)/256.0_dp !We are exploring only 50% deviation interval - and say the percentage
            !deviation is 0.5 if nothing is found
            !If we do not find a value of implie consumption there, it does not really matter if we find it
            !elsewhere.
            EPS_c = 1.0E-7_dp !tolerance level (for width of intervals).
            !ETA_c = 1.0E-12_dp !absolute value of f_x less than this accepted as zero
            ETA_c = 0.0_dp
            A_c = 0.0_dp;B_c = 0.0_dp;
            iuser = 0
            ruser = 0.0_dp

            call C05AUF(X_c, H_c, EPS_c, ETA_c, c1_imp_err, A_c, B_c, IUSER, RUSER, IFAIL)

            !testing if we really found a zero. If not, the deviation will be huge.
            if(abs(c1_imp_err(X_c,iuser,ruser)) < 1.0E-6_dp) then
                !compute the actual deviation (as percentage). Add it to agent one's residual (average)
                resid_c_all_ag(1,ok_counter) = resid_c_all_ag(1,ok_counter) + (abs(X_c - c_t1(m,1))/abs(c_t1(m,1)))/real(M_par,dp)
            else
                !otherwise the deviation is assumed to be larger than 0.5 (50 percent) because that's
                !the interval where we were looking.
                resid_c_all_ag(1,ok_counter) = resid_c_all_ag(1,ok_counter) + 0.5_dp/real(M_par,dp)
            end if

        end do

    write(*,*) 'residual at point = ',resid_all(ok_counter)
    write(*,*) 'period1 best value = ',best_t1
    write(*,*) 'period1 consumption = ',c_t1
    write(*,*) 'period1 a_prime = ',a_prime_t1
    write(*,*) 'period1 rho_prime = ',rho_prime_t1
    write(*,*) '____________________________________'

    333 end do !end of loop over sample points.


    !Save the vector of EE residuals for plotting (can do a nice historgram in Matlab)
    if(ok_counter > 0) then
    file_path = 'results/'//trim(folder_name)//'/EEresid.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) resid_all(1:ok_counter)
    close(unit = 20)
    end if
    if(ok_counter > 0) then
    file_path = 'results/'//trim(folder_name)//'/EEresid_ag1.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) resid_all_ag(1:1,1:ok_counter)
    close(unit = 20)
    end if
    if(ok_counter > 0) then
    file_path = 'results/'//trim(folder_name)//'/EEresid_ag2.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) resid_all_ag(2:2,1:ok_counter)
    close(unit = 20)
    end if

    write(string1,*) 'Share of points used in accuracy test = ',real(ok_counter,dp)/real(pars%acc_samples,dp)
    call double_output(folder_name,trim(string1))

    !compute average residual, residuals per agent, and cons interpretation of residual of agent 1
    resid_avg = 0.0_dp
    resid_avg_ag = 0.0_dp
    resid_c1_avg = 0.0_dp
    do i = 1,ok_counter
        resid_avg = resid_avg + resid_all(i) / real(ok_counter,dp)
        resid_avg_ag(1) = resid_avg_ag(1) + resid_all_ag(1,i) / real(ok_counter,dp)
        resid_avg_ag(2) = resid_avg_ag(2) + resid_all_ag(2,i) / real(ok_counter,dp)
        !Also compute the average residual of agent 1 in terms of consumption. But that is
        !not that useful because it may depend on the threshold - rather use quantiles for this purpose...
        !resid_c1_avg = resid_c1_avg + EE_c1_all(i) / real(ok_counter,dp)
        !resid_c1_avg = resid_c1_avg + resid_c_all_ag(1,ok_counter) / real(ok_counter,dp)
    end do
    resid_c1_avg = sum(reshape(resid_c_all_ag(1:1,1:ok_counter),[ok_counter]))/real(ok_counter,dp)


    write(string1,*) 'Avg EE resid = ',resid_avg
    call double_output(folder_name,trim(string1))
    write(string1,*) 'Max EE resid = ',maxval(resid_all(1:ok_counter))
    call double_output(folder_name,trim(string1))

    write(string1,*) 'Avg EE resid (agent 1) = ',resid_avg_ag(1)
    call double_output(folder_name,trim(string1))
    write(string1,*) 'Max EE resid (agent 1) = ',maxval(reshape(resid_all_ag(1,1:ok_counter),[ok_counter]))
    call double_output(folder_name,trim(string1))

    write(string1,*) 'Avg EE resid (agent 2) = ',resid_avg_ag(2)
    call double_output(folder_name,trim(string1))
    write(string1,*) 'Max EE resid (agent 2) = ',maxval(reshape(resid_all_ag(2,1:ok_counter),[ok_counter]))
    call double_output(folder_name,trim(string1))

    write(string1,*) 'Avg resid in c (agent 1) = ',resid_c1_avg
    call double_output(folder_name,trim(string1))
    write(string1,*) 'Max EE resid in c (agent 1) = ',maxval(reshape(resid_c_all_ag(1,1:ok_counter),[ok_counter]))
    call double_output(folder_name,trim(string1))
    write(string1,*) 'Min EE resid in c (agent 1) = ',minval(reshape(resid_c_all_ag(1,1:ok_counter),[ok_counter]))
    call double_output(folder_name,trim(string1))


    call double_output(folder_name,'(This is affected by choice of upper bounds, rather use quantiles...)')

    call double_output(folder_name,'_________________________________________')


    !Also compute some quantiles of residuals

    N = ok_counter
    nq = 12 !number of quantiles
    q = [0.05_dp, 0.1_dp,0.2_dp,0.3_dp,0.4_dp,0.5_dp,0.6_dp,0.7_dp,0.8_dp,0.9_dp,0.95_dp,0.99_dp] !quantile values
    qv = 0.0_dp

    ifail = -1
    call G01AMF(N,resid_all(1:ok_counter),NQ,Q,QV,IFAIL)
    !Careful - elements of resid_all are rearranged (not sure if they are all kept too - so if we want to use resid_all further,
    !we need to create a copy first (or use a copy in the quantile subroutine).
    do i = 1,nq
        write(string1,'(F4.2)') q(i)
        write(string2,*) qv(i)
        string3 = 'Quantile '//trim(string1)//' of EE resid = '//trim(string2)
        call double_output(folder_name,trim(string3))
    end do
    call double_output(folder_name,'_________________________')


    !Also get the quantiles for residuals of individual agents. At this point we overwrite resid_all
    ifail = -1
    resid_all(1:ok_counter) = reshape(resid_all_ag(1:1,1:ok_counter),[ok_counter])
    call G01AMF(N,resid_all(1:ok_counter),NQ,Q,QV,IFAIL)
    !Careful - elements of resid_all are rearranged (not sure if they are all kept too - so if we want to use resid_all further,
    !we need to create a copy first (or use a copy in the quantile subroutine).
    do i = 1,nq
        write(string1,'(F4.2)') q(i)
        write(string2,*) qv(i)
        string3 = 'Quantile '//trim(string1)//' of EE resid (agent 1) = '//trim(string2)
        call double_output(folder_name,trim(string3))
    end do
    call double_output(folder_name,'_________________________')

    ifail = -1
    resid_all(1:ok_counter) = reshape(resid_all_ag(2:2,1:ok_counter),[ok_counter])
    call G01AMF(N,resid_all(1:ok_counter),NQ,Q,QV,IFAIL)
    !Careful - elements of resid_all are rearranged (not sure if they are all kept too - so if we want to use resid_all further,
    !we need to create a copy first (or use a copy in the quantile subroutine).
    do i = 1,nq
        write(string1,'(F4.2)') q(i)
        write(string2,*) qv(i)
        string3 = 'Quantile '//trim(string1)//' of EE resid (agent 2) = '//trim(string2)
        call double_output(folder_name,trim(string3))
    end do
    call double_output(folder_name,'_________________________')

    !Also get the quantiles for residuals of individual agents. At this point we overwrite resid_all
    ifail = -1
    resid_all(1:ok_counter) = reshape(resid_c_all_ag(1:1,1:ok_counter),[ok_counter])
    call G01AMF(N,resid_all(1:ok_counter),NQ,Q,QV,IFAIL)
    !Careful - elements of resid_all are rearranged (not sure if they are all kept too - so if we want to use resid_all further,
    !we need to create a copy first (or use a copy in the quantile subroutine).
    do i = 1,nq
        write(string1,'(F4.2)') q(i)
        write(string2,*) qv(i)
        string3 = 'Quantile '//trim(string1)//' of EE resid as cons (agent 1) = '//trim(string2)
        call double_output(folder_name,trim(string3))
    end do
    call double_output(folder_name,'_________________________')



    deallocate(resid_all)
    deallocate(resid_all_ag)
    deallocate(resid_c_all_ag)
    

call double_output(folder_name,'_________________________________________')

end subroutine acc_test

!Subroutine fix_c is used to adjust a consumption vector such that it satisfies the minimum
!requirement that consumption be non-negative for all agents and all states. This is
!used as a last resort only when other methods (such as using an initial guess from previous
!period if we have confidence that this guess is good, or using multiple random initializations) failed.
!The resulting guess may be far from optimal choice and may violate other constraints - but that's a lesser
!evil than working with a consumption vector with a negative element (using logs on these leads to errors)
!or using an arbitrary constant guess.
subroutine fix_c(c,rho,pars)
    real(dp), dimension(M_par,I_par), intent(inout) :: c
    real(dp), dimension(I_par-1) :: rho !value of state rho used to set the ratio of MUs.
    type(par), intent(in) :: pars

    real(dp) :: avg_c1
    integer :: s_ind

    if(I_par > 2) then
        write(*,*) 'Error: subroutine fix_c needs to be generalized for I>2.'
    end if
    !compute the average consumption of agent 1 for all states

    avg_c1 = 0.0_dp
    do s_ind = 1,M_par
        avg_c1 = avg_c1 + c(s_ind,1)
    end do
    avg_c1 = avg_c1 / real(M_par,dp)

    !set agent 1's consumption equal across all states (to the average value)
    if(avg_c1 > 0.0_dp) then
        c(:,1) = avg_c1
    else
        c(:,1) = 0.01_dp !Consumption should never be negative on input to this subroutine but just in case...
    end if


    !And set agent 2's to rho multiple of this for all states, so that the constraint is
    !automatically satisfied.
    c(:,2) = avg_c1 * rho(1)

end subroutine fix_c

!Subroutine get_lambda computes Lagrange multipliers lambda (one per each player)
!in the social planner's problem. This is used in accuracy tests only.
!The inputs are states and functions of allocation.
subroutine get_lambda(lambda,pars,a_st,theta_lvl,g_st_ind,c,l,a_prime,R,sing_fail,omega,xi,lambda_e)
use nag_library, Only: f01blf, nag_wp
    type(par), intent(in) :: pars
    real(dp), dimension(I_par), intent(in) :: a_st
    real(dp), dimension(1,I_par), intent(in) :: theta_lvl
    !Note that we don't need rho_st
    integer :: g_st_ind !index of last-period shock realization
    real(dp), dimension(M_par,I_par), intent(in) :: c,l,a_prime
    !current consumption plan, labour plan, MU-adjusted asset holdings plan (see definition in the paper)
    real(dp), dimension(M_par,1), intent(in) :: R !interest rate plan
    real(dp), dimension(M_par,I_par), intent(out) :: lambda

    real(dp), dimension(M_par),intent(out) :: omega,xi !the multipliers omega,xi are not used to get EE residuals,
    !but they are useful for computing consumption interpretation of EE residuals. Same goes for lambda_e which is the expected
    !multiplier implied by the allocation.
    real(dp), dimension(M_par,I_par), intent(out) :: lambda_e

    logical, intent(out) :: sing_fail !on exit this is .true. if in singular value decomposition of the
    !system equation we encountered a singular value very low (which inicated near singularity). In this case we want to skip the point

    !___what follows are local variables____________
    integer :: s_ind !index for cycling over current-period shock realizations.
    real(dp), dimension(M_par,I_par) :: u_c !marginal utility of consumption
    real(dp), dimension(1,I_par) :: Eu_c !expected utility conditional on last-period shock realization

    real(dp), dimension(6,6) :: K_sys,K_sys_copy !matrix of coefficients
    real(dp), dimension(6,1) :: b_sys !matrix of RHS
    real(dp), dimension(6,1) :: x_sys !matrix of unknowns (lambda1, lambda1_e, lambda2, lambda2_e, omega, xi)

    !variables for SVD
    integer :: info,lwork
    real(dp), dimension(:), allocatable :: work

    real(dp), dimension(6) :: sigma_vec !singular values as vector (returned by NAG subroutine)
    real(dp), dimension(6,6) :: sigma_mat_inv, sigma_mat !and in matrix form (already inverted)
    real(dp), dimension(6,6) :: u_mat,VT_mat !outputs of the singular value decomposition

    logical, parameter :: debug_SVD = .false. !if this is true, the singular value decomposition accuracy
    !is tested by direct multiplication and compared to the original matrix.

    integer :: i


    sing_fail = .false.

    !First pre-compute some useful expressions
    call util_c(c,u_c,pars%A_par)
    !expected marginal utility of consumption conditional on last-period shock
    Eu_c = matmul(pars%P(g_st_ind:g_st_ind,:),u_c) / pars%A_par !devide by A_par because what we want is
    !expectation of 1/c, not A*(1/c)

    !For a current shock realization we have a system of 6 equations in 6 unknowns. We solve these
    !systems separately for each realization, and then recover the solution (lambda1, lambda2), and save it in lambda.

    !This function was done for the utility function in paper - it would need to be generalized for a general
    !utility function...

    !cycle over current shock realizations
    do s_ind = 1,M_par

    !equation(1)
    b_sys(1,1) = - pars%alpha(1,1) * pars%mass(1,1) / c(s_ind,1)

    K_sys(1,1) = pars%B_par *(l(s_ind,1) ** (1.0_dp + pars%gamma_par)) + &
        a_st(1)* ((1/(c(s_ind,1)*Eu_c(1,1)))**2)*pars%P(g_st_ind,s_ind) - 1.0_dp
    K_sys(1,2) = - (pars%Beta**2) * R(s_ind,1) * a_prime(s_ind,1)
    K_sys(1,3) = 0.0_dp
    K_sys(1,4) = 0.0_dp
    K_sys(1,5) = - (1.0_dp/theta_lvl(1,1)) * pars%B_par * l(s_ind,1) ** (pars%gamma_par)
    K_sys(1,6) = -pars%mass(1,1)

    !equation(2)
    b_sys(2,1) = - pars%alpha(1,2) * pars%mass(1,2) / c(s_ind,2)
    K_sys(2,1) = 0.0_dp
    K_sys(2,2) = 0.0_dp
    K_sys(2,3) = pars%B_par *(l(s_ind,2) ** (1.0_dp + pars%gamma_par)) + &
        a_st(2)* ((1/(c(s_ind,2)*Eu_c(1,2)))**2)*pars%P(g_st_ind,s_ind) - 1.0_dp
    K_sys(2,4) = - (pars%Beta**2) * R(s_ind,1) * a_prime(s_ind,2)
    K_sys(2,5) =  (1.0_dp/theta_lvl(1,2)) * pars%B_par * l(s_ind,2) ** (pars%gamma_par)
    K_sys(2,6) = -pars%mass(1,2)

    !equation(3)
    b_sys(3,1) = pars%alpha(1,1) * pars%mass(1,1) * pars%B_par * (l(s_ind,1) ** (pars%gamma_par))

    K_sys(3,1) = pars%B_par * c(s_ind,1) * (1.0_dp + pars%gamma_par) * (l(s_ind,1) ** (pars%gamma_par))
    K_sys(3,2) = 0.0_dp
    K_sys(3,3) = 0.0_dp
    K_sys(3,4) = 0.0_dp
    K_sys(3,5) = - (1.0_dp/theta_lvl(1,1)) * pars%B_par * (l(s_ind,1) ** (pars%gamma_par - 1.0_dp)) &
        * pars%gamma_par * c(s_ind,1)
    K_sys(3,6) = pars%mass(1,1) * theta_lvl(1,1)

    !equation(4)
    b_sys(4,1) = pars%alpha(1,2) * pars%mass(1,2) * pars%B_par * (l(s_ind,2) ** (pars%gamma_par))

    K_sys(4,1) = 0.0_dp
    K_sys(4,2) = 0.0_dp
    K_sys(4,3) = pars%B_par * c(s_ind,2) * (1.0_dp + pars%gamma_par) * (l(s_ind,2) ** (pars%gamma_par))
    K_sys(4,4) = 0.0_dp
    K_sys(4,5) = (1.0_dp/theta_lvl(1,2)) * pars%B_par * (l(s_ind,2) ** (pars%gamma_par - 1.0_dp)) &
        * pars%gamma_par * c(s_ind,2)
    K_sys(4,6) = pars%mass(1,2) * theta_lvl(1,2)

    !equation(5)
    b_sys(5,1) = 0.0_dp

    K_sys(5,1) = 1.0_dp
    K_sys(5,2) = - pars%beta * R(s_ind,1)
    K_sys(5,3) = 0.0_dp
    K_sys(5,4) = 0.0_dp
    K_sys(5,5) = 0.0_dp
    K_sys(5,6) = 0.0_dp

    !equation(6)
    b_sys(6,1) = 0.0_dp
    K_sys(6,1) = 0.0_dp
    K_sys(6,2) = 0.0_dp
    K_sys(6,3) = 1.0_dp
    K_sys(6,4) = - pars%beta * R(s_ind,1)
    K_sys(6,5) = 0.0_dp
    K_sys(6,6) = 0.0_dp


    !Solve the system by matrix inversion - it should be fairly non-singular so no issues here
    !But also check and display a warning if singularity is an issue.

    !K_sys = 0.0_dp
    do i=1,6
        K_sys(i,i) = K_sys(i,i) + 0.0_dp
    end do

    !copy K_sys (it will be overwritten by its inverse)
    K_sys_copy = K_sys

    !Using singular value decomposition to solve the system (and detect problems)
    lwork = max(1,3*min(6,6) + max(6,6),5*min(6,6))
    lwork = lwork + min(6,6)*1024 !the constant should be optimal block size.
    !(but in the case of small matrix such as here, we could just use a large constant - this code
    !was taken from TaxHetGrReg where it's used in a more complex setting).
    allocate(work(lwork))

    call F08KBF('S','S',6,6,K_sys,6,sigma_vec,U_mat,6,VT_mat,6,WORK,LWORK,INFO)
    !Careful - matrix VT_mat already contains the transposed matrix V (see documentation of the subroutine)

    !Matrix K_sys was changed by the subroutine for SVD! After this point we can only use K_sys_copy
    deallocate(work)

    if(info>0.0_dp) then
        write(*,*) 'Error: SVD in subroutine reg_M_prep failed!'
        error stop
    end if

    !Get the matrix form of sigma (and its inverse)
    sigma_mat = 0.0_dp
    sigma_mat_inv = 0.0_dp
    do i=1,6
        sigma_mat(i,i) = sigma_vec(i) !actually we need
        !this only if SVD_debug == .true., but it's a small object so it doesn't matter that we get it
        !all the time.
        sigma_mat_inv(i,i) = 1.0_dp/sigma_vec(i)
    end do

    !If some singular value is very small (less than 0.00001) then the system is pretty close to singular.
    !display a warning (if this happens often, then we want to introduce another argument of this subroutine),
    !so that we can find out outside of the subroutine whether this is the case, and discard the point in state
    !space at which this happens from the accuracy test.
    if(minval(abs(sigma_vec)) < 0.00001_dp) then
        !write(*,*) 'Singular value is (in absolute value) minimally ',minval(abs(sigma_vec))
        sing_fail = .true.
        return
    end if

    if(DEBUG_SVD) then
        write(*,*) 'Debug_SVD in subroutine get_lambda (hardcoded parameter).'
        write(*,*) 'If SVD done properly, the following two numbers should be pretty much exactly zero.'
        !compute the multiplication which should yield the original matrix K (if SVD is correct),
        !then subtract K_sys, and write it out on the screen - should be zeroes (pretty much exactly)
        write(*,*) 'Error in SVD = '
        write(*,*) maxval(reshape(abs(K_sys_copy - matmul(U_mat, matmul(sigma_mat,VT_mat))),[36]))
        write(*,*) 'Error in getting inverse by SVD and multiplying the original matrix = '

        !Get the inverse of matrix K using the SVD, and save it in K_sys.
        K_sys = matmul(transpose(VT_mat),matmul(sigma_mat_inv,transpose(U_mat)))
        K_sys = matmul(K_sys,K_sys_copy) !Multiply by the original matrix, should get unit matrix, save it in K_sys
        !subtract the ones
        do i=1,6
            K_sys(i,i) = K_sys(i,i) - 1.0_dp
        end do
        !the result should be close to zero (max abs value)
        write(*,*) maxval(reshape(abs(K_sys),[36]))
        error stop
    end if

    !Obtain solution of the system (fairly obvious linear algebra - check subroutine for SVD in case of doubt)
    !(as column vector)
    x_sys = matmul(matmul(transpose(VT_mat),matmul(sigma_mat_inv,transpose(U_mat))),b_sys)

    lambda(s_ind,1) = x_sys(1,1) !lambda of agent 1 in state s_ind
    lambda(s_ind,2) = x_sys(3,1) !lambda of agent 2 in state s_ind

    !The following multipliers are not used in computing EE residuals, but are used for getting an
    !interpretation of EE residuals in terms of consumption.
    lambda_e(s_ind,1) = x_sys(2,1)
    lambda_e(s_ind,2) = x_sys(4,1)
    omega(s_ind) = x_sys(5,1)
    xi(s_ind) = x_sys(6,1)


    end do



end subroutine get_lambda

!Function c1_imp_err computes error as function of implied consumption c1_imp. This is error in the non-linear equation
!in c1(st), and if this error is zero, the equation is solved. If we want to compute this error for agent 2, then we
!have to change this slightly. c1_imp and lambda_t1_imp are scalars, corresponding to one state only.
!
!Warning: this function already presumes the logarithmic utility function. We will need to generalize this if
!we are interested in accuracy tests for other functions.
function c1_imp_err(x,iuser,ruser)
    real(dp) :: c1_imp_err
    integer :: iuser(*)
    real(dp) :: x,ruser(*)


    real(dp), dimension(M_par,I_par) :: c_pl_a !consumption plan (for all agents at all states)
    integer :: s_ind_a !index of current-period shock realization.
    integer :: s_st_ind_a !index of last-period shock realization.
    real(dp) :: A_a,B_a,C_a !three pre-computed terms
    real(dp), dimension(M_par,M_par) :: P_a !probability transition matrix

    real(dp) :: Ec1recip !expected reciprocal of agent 1's consumption.
    integer :: m


    common /c1_imp_err_comm/ c_pl_a,A_a,B_a,C_a,P_a,s_st_ind_a,s_ind_a

    c1_imp_err = 0.0_dp !initialize value

    !in c_pl, change the element (corresponding to agent 1's consumption in state s_ind) to the implied consumption
    !for which we are trying to compute the error
    c_pl_a(s_ind_a,1) = x

    !Ec1recip
    Ec1recip = 0.0_dp
    do m = 1,M_par
        Ec1recip = Ec1recip + (1.0_dp/c_pl_a(m,1)) * P_a(s_st_ind_a,m)
    end do

    !The error in equation (remember: x = c_pl_a(s_ind_a,1))
    c1_imp_err = A_a/x + B_a * ((x*Ec1recip)**(-2.0_dp)) + C_a

end function c1_imp_err


end module mod_taxhetgr

