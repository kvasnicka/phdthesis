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

!type grids_type contains grids. It is written for 2 agents only because
!the addition of third player would make the problem computationally infeasible
!using currently available hardware (going from 5 state variables to 7)
type grids_type
    real(dp),allocatable,dimension(:) :: A1_gr !asset holdings of agent 1
    !real(dp),allocatable,dimension(:) :: A2_gr !asset holdings of agent 2
    real(dp),allocatable,dimension(:) :: rho_gr !u1/u2
    real(dp),allocatable,dimension(:) :: T_gr !relative productivity
end type grids_type

!type FP_result_type contains all sorts of results of the first-period problem
type FP_result_type
    real(dp), dimension(1,I_par) :: c,l,u_c,u_l
    real(dp), dimension(1,I_par-1) :: b_net_prime,a_prime
    real(dp) :: rho_prime !will need to increase dimension for I>2
    real(dp) :: t_prime !'time variable'
    integer :: s_ind_init !initial index - this is usually drawn from a stationary distribution
    !(so that we don't get a solution dependent on initial state which can imply large quantitative
    !trend in government expenditure in expectation).

end type FP_result_type

!Type sim_result type contains
type sim_result_type
    !First index is time period, 2nd is agent
    real(dp), allocatable, dimension(:,:) :: c_sim,l_sim,b_sim,b_prime_sim,b_net_sim,b_net_prime_sim,a_prime_sim,u_c_sim,u_l_sim,theta_lvl_sim
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
    real(dp), allocatable, dimension(:) :: Y_sim,B_gov_sim,B_gov_to_GDP_sim,trans_sim,trans_to_GDP_sim
end type sim_result_type


!testing (maybe it will not work properly with coarrays)
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

    real(dp) :: tau_first,tau_last,tau_lastminfirst 


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

    !and save quantiles of tax policy
    real(dp), dimension(:), allocatable :: tau_05, tau_95
    real(dp), dimension(:), allocatable :: trans_05, trans_95


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

    call sim_result_alloc(sim_plus_sim,pars_comm,1)


    sim_plus_sim%c_sim = sr1%c_sim + sr2%c_sim
    sim_plus_sim%l_sim = sr1%l_sim + sr2%l_sim
    sim_plus_sim%b_sim = sr1%b_sim + sr2%b_sim
    sim_plus_sim%b_prime_sim = sr1%b_prime_sim + sr2%b_prime_sim
    sim_plus_sim%b_net_sim = sr1%b_net_sim + sr2%b_net_sim
    sim_plus_sim%b_net_prime_sim = sr1%b_net_prime_sim + sr2%b_net_prime_sim
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
    sim_plus_sim%trans_sim = sr1%trans_sim + sr2%trans_sim

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
    sim_min_sim%b_net_sim = sr1%b_net_sim - sr2%b_net_sim
    sim_min_sim%b_net_prime_sim = sr1%b_net_prime_sim - sr2%b_net_prime_sim
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
    sim_min_sim%trans_sim = sr1%trans_sim - sr2%trans_sim

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
    sim_multi%b_net_sim = sr%b_net_sim * multi
    sim_multi%b_net_prime_sim = sr%b_net_prime_sim * multi
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
    sim_multi%trans_sim = sr%trans_sim * multi
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
            sim_result%b_net_sim(pars%T_sim,I_par-1),sim_result%b_net_prime_sim(pars%T_sim,I_par-1),sim_result%a_prime_sim(pars%T_sim,I_par-1),&
            sim_result%u_c_sim(pars%T_sim,I_par),sim_result%b_prime_sim(pars%T_sim,I_par),&
            sim_result%u_l_sim(pars%T_sim,I_par),sim_result%theta_lvl_sim(pars%T_sim,I_par))
        allocate(sim_result%rho_prime_sim(pars%T_sim),sim_result%t_prime_sim(pars%T_sim),&
            sim_result%t_sim(pars%T_sim),sim_result%s_ind_sim(pars%T_sim),sim_result%s_real_sim(pars%T_sim),&
            sim_result%R_sim(pars%T_sim-1),sim_result%tau_sim(pars%T_sim),sim_result%g_lvl_sim(pars%T_sim),&
            sim_result%gini_income_sim(pars%T_sim),sim_result%gini_prod_sim(pars%T_sim))
        allocate(sim_result%Y_sim(pars%T_sim),sim_result%B_gov_sim(pars%T_sim),&
            sim_result%B_gov_to_GDP_sim(pars%T_sim),sim_result%trans_to_GDP_sim(pars%T_sim))
        allocate(sim_result%G_to_Y_sim(pars%T_sim),sim_result%gini_allincome_sim(pars%T_sim),sim_result%trans_sim(pars%T_sim))

    else if(wtd==0) then
        !Deallocate
        deallocate(sim_result%c_sim,sim_result%l_sim,sim_result%b_sim,sim_result%b_net_sim,sim_result%b_net_prime_sim,sim_result%a_prime_sim,&
            sim_result%u_c_sim,sim_result%u_l_sim,sim_result%theta_lvl_sim,sim_result%b_prime_sim)
        deallocate(sim_result%rho_prime_sim,sim_result%t_prime_sim,sim_result%t_sim,&
            sim_result%s_ind_sim,sim_result%s_real_sim,sim_result%R_sim,sim_result%tau_sim,sim_result%g_lvl_sim,&
            sim_result%gini_income_sim,sim_result%gini_prod_sim)
        deallocate(sim_result%Y_sim,sim_result%B_gov_sim,sim_result%B_gov_to_GDP_sim,sim_result%trans_to_GDP_sim)
        deallocate(sim_result%G_to_Y_sim,sim_result%gini_allincome_sim,sim_result%trans_sim)
    else
        write(*,*) 'Error: wrong value of wtd in subroutine sim_result_alloc.'
        sync all
        error stop
    end if


end subroutine sim_result_alloc


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

!Subroutine check_constr checks whether constraints are violated. It also computes a convex
!cost of violating constraints (so far this is equal to 0)
subroutine check_constr(c,l,l_max_mat,theta_lvl,a_prime,a_min_mat,a_max_mat&
,what_binds,rho_prime,rho_min_vec,rho_max_vec,constr_ok,loss,A_par,rofv)
    real(dp), dimension(M_par,I_par), intent(in) :: c,l
    real(dp), dimension(M_par,I_par-1), intent(in) :: a_prime
    real(dp), dimension(M_par,I_par), intent(in) :: l_max_mat
    real(dp), dimension(M_par,I_par-1), intent(in) :: a_min_mat,a_max_mat
    integer, dimension(3), intent(in) :: what_binds
    real(dp), dimension(M_par,1), intent(in) :: rho_min_vec,rho_max_vec,rho_prime
    real(dp), dimension(1,I_par), intent(in) :: theta_lvl
    real(dp), intent(in) :: A_par


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


    real(dp), dimension(M_par,1) :: ones = 1.0_dp

    real(dp), dimension(M_par,I_par) :: u_c !marginal utility of consumption

    real(dp), dimension(M_par,I_par) :: b_min_mat_local, b_max_mat_local !Local copies: These will
    !be used to check constraints (will be essentially the same as the matrices passed to this subroutine
    !possibly with one constraint replaced by a large constant).

    real(dp) :: loss1a,loss2a,loss3a !levels of loss from least severe problem to the most severe


!The function checks all constraints - even if one is violated we still go through till the end
!and add a cost of violation which should make it unappealing for the program to choose such point.
!
!There will be two loss conditions: a strict one for violating a fundamental constraint, and
!a milder one for violating constrain on marginal-utility adjusted


loss1a = 10.0_dp
loss2a = 100.0_dp
loss3a = 10000.0_dp


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

!write(*,*) 'loss after checking c: =',loss

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


!Check constraints on rho
!First need to compute the marginal utility of consumption from c and l

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


end subroutine check_constr
!
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

        real(dp), dimension(:,:,:,:), pointer :: V_old_pntr !pointer to V_old, passed using common block

        !The following variables are passed to the function using a common block
        real(dp), dimension(:), pointer :: a1_gr_pntr,rho_gr_pntr,t_gr_pntr
        real(dp), dimension(:,:), pointer :: S_pntr,P_pntr, mass_pntr
        !and transition matrix.
        real(dp), pointer :: A_par_pntr,B_par_pntr,Gamma_par_pntr
        type(par), pointer :: pars_pntr

        !Values of state variables as saved in the grid (hence the suffix st).
        real(dp), dimension(I_par-1) :: a_st
        real(dp), dimension(I_par-1) :: rho_st
        real(dp) :: t_st,g_st
        integer :: g_st_ind
        real(dp), dimension(1,I_par) :: theta_lvl

        real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat
        real(dp), dimension(M_par,I_par-1) :: a_min_mat,a_max_mat
        real(dp), dimension(M_par,1) :: rho_min_vec,rho_max_vec

        real(dp), dimension(M_par,I_par) :: g_lvl_mat
        real(dp), dimension(M_par,M_par) :: P_onerow

        !Local variables
        real(dp), dimension(M_par,I_par) :: c_inc_mat = 0.0_dp !incomplete consumption matrix
        !Labour
        real(dp), dimension(M_par,I_par) :: lab_mat = 0.0_dp
        !next-period states
        real(dp) :: t_st_pr
        real(dp), dimension(M_par,I_par-1) :: a_prime!,b_prime

        logical :: c_im_found = .false., constr_ok

        integer :: g_ind !index of shock realization

        real(dp) :: loss !loss from violating constraints.

        real(dp), dimension(M_par,1) :: conval_all !continuation value in all states

        !rho_prime and u_c
        real(dp), dimension(M_par,I_par) :: u_c !marginal utility of consumption
        real(dp), dimension(M_par,1) :: rho_prime

        !temporary variable for debugging
        real(dp), dimension(M_par,M_par) :: P_debug

        common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
        !to some subroutines.
        common /pars_pointer/ pars_pntr

        common /cpc_ret_bloc/ a_st,theta_lvl,rho_st,t_st,g_st,g_st_ind
        common /cpc_ret_bloc2/ g_lvl_mat, P_onerow,l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec

        common /comm1/ a1_gr_pntr,rho_gr_pntr,t_gr_pntr,A_par_pntr,B_par_pntr,&
        Gamma_par_pntr,S_pntr,P_pntr,mass_pntr

        common /commloss/ loss

        !Get the incomplete consumption matrix from the input vector (save it as c_inc_mat)
        call c_inc_vtm(xc,c_inc_mat)


        call get_c_im(c_inc_mat,iuser,a_st,rho_st,theta_lvl,g_lvl_mat(:,1),g_st_ind,pars_pntr,P_onerow,c_im_found)
        !(from this point on, c_inc_mat is actually the complete consumption plan)

        if(.not. c_im_found) then
            fc = 10000.0_dp !significant discontinuity, hopefully will almost never happen
            !because it could mess with the minimization algorithm. Also perhaps this could
            !be less...
            return
        end if


        !Compute Labour (for that we need the state passed in a common block)
        call lab(c_inc_mat,lab_mat,theta_lvl,mass_pntr,g_lvl_mat(:,1),Gamma_par_pntr)


        !Now that we know consumption and labour plan, we can compute the current (expected) return
        !Note: will need to pass mass separately and get rid of pars_pntr.
        fc = SWF(c_inc_mat,lab_mat,g_st_ind,P_pntr,mass_pntr,pars_pntr%alpha,A_par_pntr,B_par_pntr,Gamma_par_pntr)

        !Compute next-period states so we can get continuation value.
        t_st_pr = t_prime(t_st,t_gr_pntr(pars_pntr%N_t))

        call ass(a_st,c_inc_mat,lab_mat,pars_pntr%A_par,pars_pntr%B_par,pars_pntr%Gamma_par,&
        pars_pntr%beta,P_onerow,a_prime)

        !In TaxHetGrTr we don't use constraints on b yet
        !        call get_b_prime(b_prime,a_prime,c_inc_mat,pars_pntr)

        call get_rho_prime(c_inc_mat,rho_prime)


        constr_ok = .false.
        !Check constraints (Get a loss value). All of the matrices/vectors for checking constraints need to
        !be precomputed and passed to this thing by common block.


        call check_constr(c_inc_mat,lab_mat,l_max_mat,theta_lvl,a_prime,a_min_mat&
        ,a_max_mat,iuser,rho_prime,rho_min_vec,rho_max_vec,constr_ok,loss,pars_pntr%A_par,.false.)
        !Now add discounted expected continuation value: This will be computed by function dECV


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
            call interp_V(V_old_pntr,[a_prime(g_ind,1),rho_prime(g_ind,1),t_st_pr]&
            ,g_ind,conval_all(g_ind,1),a1_gr_pntr,rho_gr_pntr,t_gr_pntr,size(a1_gr_pntr),&
            size(rho_gr_pntr),size(t_gr_pntr))
        end do

        !Now we have continuation value for every possible shock realization. Add the discounted
        !expectation of this term to the current return. (g_st_ind is the last shock realization)

        fc = fc + pars_pntr%beta*maxval(matmul(P_pntr(g_st_ind,:),conval_all))
        !maxval is there simply to convert (1,1) array to scalar. Maybe this can be done in a better way...

        !Subtract loss stemming from violating some constraints.
        fc = (fc - loss)*(-1.0_dp)
        !We need (-1) multiple of return here (minimization subroutine)!!!

end subroutine CPC_ret

!This is just a copy of CPC_ret with one less input - so it can be used by a 1-d optimization
!subroutine. In the future it might be better to have this just as an 'interface subroutine'
!which would call CPC_ret! Much better for making changes, avoiding bugs, etc.
subroutine CPC_ret_1d(xc,fc,iuser,ruser)
        real (dp), intent (out) :: fc
        real (dp), intent (inout) :: ruser(*)
        real (dp), intent (in) :: xc
        integer, intent (inout) :: iuser(*) !In this instance, it is used to pass what_binds

        real(dp), dimension(:,:,:,:), pointer :: V_old_pntr !pointer to V_old, passed using common block

        !The following variables are passed to the function using a common block
        real(dp), dimension(:), pointer :: a1_gr_pntr,rho_gr_pntr,t_gr_pntr
        real(dp), dimension(:,:), pointer :: S_pntr,P_pntr, mass_pntr
        !and transition matrix.
        real(dp), pointer :: A_par_pntr,B_par_pntr,Gamma_par_pntr
        type(par), pointer :: pars_pntr

        !Values of state variables as saved in the grid (hence the suffix st).
        real(dp), dimension(I_par-1) :: a_st
        real(dp), dimension(I_par-1) :: rho_st
        real(dp) :: t_st,g_st
        integer :: g_st_ind
        real(dp), dimension(1,I_par) :: theta_lvl

        real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat
        real(dp), dimension(M_par,I_par-1) :: a_min_mat,a_max_mat
        real(dp), dimension(M_par,1) :: rho_min_vec,rho_max_vec

        real(dp), dimension(M_par,I_par) :: g_lvl_mat
        real(dp), dimension(M_par,M_par) :: P_onerow

        !Local variables
        real(dp), dimension(M_par,I_par) :: c_inc_mat = 0.0_dp !incomplete consumption matrix
        !Labour
        real(dp), dimension(M_par,I_par) :: lab_mat = 0.0_dp
        !next-period states
        real(dp) :: t_st_pr
        real(dp), dimension(M_par,I_par-1) :: a_prime!,b_prime

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

        common /comm1/ a1_gr_pntr,rho_gr_pntr,t_gr_pntr,A_par_pntr,B_par_pntr,&
        Gamma_par_pntr,S_pntr,P_pntr,mass_pntr

        common /commloss/ loss

        !Get the incomplete consumption matrix from the input vector (save it as c_inc_mat)
        !call c_inc_vtm(xc,c_inc_mat)
        c_inc_mat = xc !pretty much the only change from CPC_ret

        call get_c_im(c_inc_mat,iuser,a_st,rho_st,theta_lvl,g_lvl_mat(:,1),g_st_ind,pars_pntr,P_onerow,c_im_found)
        !(from this point on, c_inc_mat is actually the complete consumption plan)

        if(.not. c_im_found) then
            fc = 10000.0_dp !significant discontinuity, hopefully will almost never happen
            !because it could mess with the minimization algorithm. Also perhaps this could
            !be less...
            return
        end if


        !Compute Labour (for that we need the state passed in a common block)
        call lab(c_inc_mat,lab_mat,theta_lvl,mass_pntr,g_lvl_mat(:,1),Gamma_par_pntr)


        !Now that we know consumption and labour plan, we can compute the current (expected) return
        !Note: will need to pass mass separately and get rid of pars_pntr.
        fc = SWF(c_inc_mat,lab_mat,g_st_ind,P_pntr,mass_pntr,pars_pntr%alpha,A_par_pntr,B_par_pntr,Gamma_par_pntr)


        !Compute next-period states so we can get continuation value.
        t_st_pr = t_prime(t_st,t_gr_pntr(pars_pntr%N_t))

        call ass(a_st,c_inc_mat,lab_mat,pars_pntr%A_par,pars_pntr%B_par,pars_pntr%Gamma_par,&
        pars_pntr%beta,P_onerow,a_prime)

        !In TaxHetGrTr we don't use constraints on b yet
        !        call get_b_prime(b_prime,a_prime,c_inc_mat,pars_pntr)

        call get_rho_prime(c_inc_mat,rho_prime)


        constr_ok = .false.
        !Check constraints (Get a loss value). All of the matrices/vectors for checking constraints need to
        !be precomputed and passed to this thing by common block.

        call check_constr(c_inc_mat,lab_mat,l_max_mat,theta_lvl,a_prime,a_min_mat&
        ,a_max_mat,iuser,rho_prime,rho_min_vec,rho_max_vec,constr_ok,loss,pars_pntr%A_par,.false.)
        !Now add discounted expected continuation value: This will be computed by function dECV


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
            call interp_V(V_old_pntr,[a_prime(g_ind,1),rho_prime(g_ind,1),t_st_pr]&
            ,g_ind,conval_all(g_ind,1),a1_gr_pntr,rho_gr_pntr,t_gr_pntr,size(a1_gr_pntr),&
            size(rho_gr_pntr),size(t_gr_pntr))
        end do

        !Now we have continuation value for every possible shock realization. Add the discounted
        !expectation of this term to the current return. (g_st_ind is the last shock realization)

        fc = fc + pars_pntr%beta*maxval(matmul(P_pntr(g_st_ind,:),conval_all))
        !maxval is there simply to convert (1,1) array to scalar. Maybe this can be done in a better way...

        !Subtract loss stemming from violating some constraints.
        fc = (fc - loss)*(-1.0_dp)
        !We need (-1) multiple of return here (minimization subroutine)!!!


end subroutine CPC_ret_1d


!!Subroutine CPC_ret2 is a copy of CPC_ret, the only difference being the presence of one
!!more output. This is defined for use in a different NAg subroutine than CPC_ret.

subroutine CPC_ret2(n,xc,fc,iuser,ruser,inform)
        real (dp), intent (out) :: fc
        integer, intent (in) :: n
        real (dp), intent (inout) :: ruser(*)
        real (dp), intent (in) :: xc(n)
        integer, intent (inout) :: iuser(*) !In this instance, it is used to pass what_binds
        integer, intent(out) :: inform

        real(dp), dimension(:,:,:,:), pointer :: V_old_pntr !pointer to V_old, passed using common block

        !The following variables are passed to the function using a common block
        real(dp), dimension(:), pointer :: a1_gr_pntr,rho_gr_pntr,t_gr_pntr
        real(dp), dimension(:,:), pointer :: S_pntr,P_pntr, mass_pntr
        !and transition matrix.
        real(dp), pointer :: A_par_pntr,B_par_pntr,Gamma_par_pntr
        type(par), pointer :: pars_pntr

        !Values of state variables as saved in the grid (hence the suffix st).
        real(dp), dimension(I_par-1) :: a_st
        real(dp), dimension(I_par-1) :: rho_st
        real(dp) :: t_st,g_st
        integer :: g_st_ind
        real(dp), dimension(1,I_par) :: theta_lvl

        real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat
        real(dp), dimension(M_par,I_par-1) :: a_min_mat,a_max_mat
        real(dp), dimension(M_par,1) :: rho_min_vec,rho_max_vec

        real(dp), dimension(M_par,I_par) :: g_lvl_mat
        real(dp), dimension(M_par,M_par) :: P_onerow

        !Local variables
        real(dp), dimension(M_par,I_par) :: c_inc_mat = 0.0_dp !incomplete consumption matrix
        !Labour
        real(dp), dimension(M_par,I_par) :: lab_mat = 0.0_dp
        !next-period states
        real(dp) :: t_st_pr
        real(dp), dimension(M_par,I_par-1) :: a_prime!,b_prime

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

        common /comm1/ a1_gr_pntr,rho_gr_pntr,t_gr_pntr,A_par_pntr,B_par_pntr,&
        Gamma_par_pntr,S_pntr,P_pntr,mass_pntr

        common /commloss/ loss

        inform = 0 !This is very important because its absence leads to bugs when the number of grid points
        !is relatively large (memory access issues)

        !Get the incomplete consumption matrix from the input vector (save it as c_inc_mat)
        call c_inc_vtm(xc,c_inc_mat)


        call get_c_im(c_inc_mat,iuser,a_st,rho_st,theta_lvl,g_lvl_mat(:,1),g_st_ind,pars_pntr,P_onerow,c_im_found)
        !(from this point on, c_inc_mat is actually the complete consumption plan)

        if(.not. c_im_found) then
            fc = 10000.0_dp !significant discontinuity, hopefully will almost never happen
            !because it could mess with the minimization algorithm. Also perhaps this could
            !be less...
            return
        end if


        !Compute Labour (for that we need the state passed in a common block)
        call lab(c_inc_mat,lab_mat,theta_lvl,mass_pntr,g_lvl_mat(:,1),Gamma_par_pntr)


        !Now that we know consumption and labour plan, we can compute the current (expected) return
        !Note: will need to pass mass separately and get rid of pars_pntr.
        fc = SWF(c_inc_mat,lab_mat,g_st_ind,P_pntr,mass_pntr,pars_pntr%alpha,A_par_pntr,B_par_pntr,Gamma_par_pntr)


        !Compute next-period states so we can get continuation value.
        t_st_pr = t_prime(t_st,t_gr_pntr(pars_pntr%N_t))

        call ass(a_st,c_inc_mat,lab_mat,pars_pntr%A_par,pars_pntr%B_par,pars_pntr%Gamma_par,&
        pars_pntr%beta,P_onerow,a_prime)

        !In TaxHetGrTr we don't use constraints on b yet
        !        call get_b_prime(b_prime,a_prime,c_inc_mat,pars_pntr)

        call get_rho_prime(c_inc_mat,rho_prime)


        constr_ok = .false.
        !Check constraints (Get a loss value). All of the matrices/vectors for checking constraints need to
        !be precomputed and passed to this thing by common block.

        call check_constr(c_inc_mat,lab_mat,l_max_mat,theta_lvl,a_prime,a_min_mat&
        ,a_max_mat,iuser,rho_prime,rho_min_vec,rho_max_vec,constr_ok,loss,pars_pntr%A_par,.false.)
        !Now add discounted expected continuation value: This will be computed by function dECV

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
            call interp_V(V_old_pntr,[a_prime(g_ind,1),rho_prime(g_ind,1),t_st_pr]&
            ,g_ind,conval_all(g_ind,1),a1_gr_pntr,rho_gr_pntr,t_gr_pntr,size(a1_gr_pntr),&
            size(rho_gr_pntr),size(t_gr_pntr))
        end do

        !Now we have continuation value for every possible shock realization. Add the discounted
        !expectation of this term to the current return. (g_st_ind is the last shock realization)

        fc = fc + pars_pntr%beta*maxval(matmul(P_pntr(g_st_ind,:),conval_all))
        !maxval is there simply to convert (1,1) array to scalar. Maybe this can be done in a better way...

        !Subtract loss stemming from violating some constraints.
        fc = (fc - loss)*(-1.0_dp)
        !We need (-1) multiple of return here (minimization subroutine)!!!


end subroutine CPC_ret2


!Subroutine CTS_num_points computes the number of grid points at an iteration of CTS algorithm.
subroutine CTS_num_points(CTS_ind,pars,N_a,N_rho,N_T)
    integer, intent(in) :: CTS_ind
    type(par), intent(in) :: pars
    integer, intent(out) :: N_a,N_rho,N_T

    integer :: split = 0 !times that the number of gridpoints is split at this iteration of CTS algorithm

    if(pars%CTS_split_times>0) then
        split = pars%CTS_split_times + 1 - CTS_ind

        N_a = floor(pars%N_a * 0.5_dp**split);
        N_rho = floor(pars%N_rho * 0.5_dp**split);
        N_T = floor(pars%N_T * 0.5_dp**split);

        !If some resulting number of grid points is lower than the minimum number of grid points
        !per that variable, adjust this. However, do this only if pars%CTS_split_times > 0. If this is
        !zero it means that we are not using CTS algorithm and it makes no sense to use the floor.
        N_a=maxval([pars%CTS_gridpoints_floor,N_a])
        N_rho=maxval([pars%CTS_gridpoints_floor,N_rho])
        N_T=maxval([pars%CTS_gridpoints_floor,N_T])
    else !CTS algorithm not used, use the number of gridpoints given in parameter file
        N_a = pars%N_a
        N_rho = pars%N_rho
        N_T = pars%N_T
    end if

    !If static_environment is used, it means that there is no growth. In this case
    !set N_theta = 1.
    if(pars%static_environment) N_T = 1

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
!
!The product of this is not fully used in the current version of the program
!where borrowing constraints were not implemented yet.
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

!This subroutine will compute b_net_prime in all states. It is used only in simulation to recover
!net asset holdings from the marginally adjusted ones that serve as state variable in V
subroutine get_b_prime(b_net_prime,a_prime,c,pars)
   real(dp), dimension(M_par,I_par-1), intent(inout) :: b_net_prime
   real(dp), dimension(M_par,I_par-1), intent(in) :: a_prime
   real(dp), dimension(M_par,I_par), intent(in) :: c
   type(par), intent(in) :: pars

   real(dp), dimension(M_par,I_par) :: u_c
   !get marginal utility of consumption
   call util_c(c,u_c,pars%A_par)

   !remember: a_i =(u_i*b_i)/beta
   !and the first column of b_net_prime and a_prime belongs to agent 2 (different convention
   !than for c and u_c, where column i belongs to agent i
   b_net_prime(:,1) = pars%beta*a_prime(:,1)/u_c(:,2)

   if(I_par>2) then
        write(*,*) 'Error: subroutine get_b_prime must be generalized for I>2.'
        error stop
   end if

end subroutine get_b_prime


!subroutine get_c_guess_img computes an initial guess of consumption for all agents in all states
!at every gridpoint on which an image works. It also computes the share of nodes at which this is feasible.
!This result is saved in coarray share_feas_img.
!The feasible guesses are returned in vector of type c_guess_all (see its definition in this module).
subroutine get_c_guess_img(c_guess_img,min_ind_img,max_ind_img,ind_unfold_all,pars,grids,share_feas_img,avg_loss_img,LFFC_choices,num_fc)
    real(dp), codimension[*] :: share_feas_img !The weighted share of points where a feasible c was found
    real(dp) :: avg_loss_img
    integer, intent(in) :: min_ind_img,max_ind_img
    type(par), codimension[*] :: pars !parameters. Not sure why this is a coarray here... Could be a coarray
    !only where we need to share data, here it could be just a standard type. Low priority issue though.
    type(c_guess_all), dimension(min_ind_img:max_ind_img), intent(out) :: c_guess_img
    integer, dimension(4,min_ind_img:max_ind_img) :: ind_unfold_all !4 states if I=2
    real(dp), dimension(:,:,:), intent(in) :: LFFC_choices
    integer, intent(in) :: num_fc
    type(grids_type), intent(in) :: grids

    integer :: gp_ind !index for cycling over gridpoints. The main index from which we will get the indices
    !for the state variables from ind_unfold_all
    !Some of the following variables will become row vectors.
    !integer :: a1_ind,a2_ind,rho_ind,theta_ind,m_ind !- indices not needed to be kept separately?

    !Values of state variables as saved in the grid (hence the suffix st).
    real(dp), dimension(I_par-1) :: a_st
    real(dp), dimension(I_par-1) :: rho_st
    real(dp) :: t_st,g_st
    integer :: g_st_ind

    !Some variables used in calling subroutine LFFC
    real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat
    real(dp), dimension(M_par,I_par-1) :: a_min_mat,a_max_mat
    real(dp), dimension(M_par,1) :: rho_min_vec,rho_max_vec
    type(c_guess_all) :: c_guess_all_gp !output: all feasible choices at a grid point
    logical :: fc_found !true if a feasible choice found at this gridpoint
    !possibilities of what_binds
    integer, dimension(1+I_par*M_par,3) :: WB_poss

    real(dp) :: rel_prod_min !minimum relative productivity. This will be the lowest relative productivity
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
        !a_st(2) = grids%a2_gr(ind_unfold_all(2,gp_ind))
        rho_st = grids%rho_gr(ind_unfold_all(2,gp_ind))
        t_st = grids%t_gr(ind_unfold_all(3,gp_ind))
        g_st = pars%S(1,ind_unfold_all(4,gp_ind))
        g_st_ind = ind_unfold_all(4,gp_ind)



        call LFFC(a_st,rho_st,t_st,g_st,g_st_ind,pars,b_min_mat,b_max_mat,l_max_mat,a_min_mat,a_max_mat,&
        rho_min_vec,rho_max_vec,LFFC_choices,num_fc,WB_poss,c_guess_all_gp,fc_found,loss_gp)




        !save the guess
        c_guess_img(gp_ind) = c_guess_all_gp
        !update statistics on feasibility of shares, etc.
        avg_loss_img = avg_loss_img + loss_gp/real(max_ind_img-min_ind_img+1)

        if(fc_found) then
            share_feas_img = share_feas_img + 1.0_dp/real(max_ind_img-min_ind_img+1)
        end if


    end do


end subroutine get_c_guess_img



!subroutine get_c_im computes element I,M of the consumption matrix, depending on what_binds,
!i.e., on which borrowing constraint binds (if any).
subroutine get_c_im(cons_mat,what_binds,a_st,rho_st,theta_lvl,g_lvl_vec,g_ind_st,pars,P_onerow,c_im_found)
    real(dp), dimension(M_par,I_par), intent(inout) :: cons_mat
    integer, dimension(3), intent(in) :: what_binds
    real(dp), dimension(I_par-1), intent(in) :: a_st
    real(dp), dimension(I_par-1), intent(in) :: rho_st
    real(dp), dimension(1,I_par), intent(in) :: theta_lvl
    real(dp), dimension(M_par,1), intent(in) :: g_lvl_vec !gov expenditure in levels in all states
    integer, intent(in) :: g_ind_st !index
    type(par), intent(in) :: pars
    real(dp), dimension(M_par,M_par), intent(in) :: P_onerow
    logical, intent(inout) :: c_im_found

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

        if(c_im_tmp>0.0_dp) then
            cons_mat(M_par,I_par) = c_im_tmp
            c_im_found = .true.
            return
        else
            cons_mat(M_par,I_par) = 0.00001_dp !just to avoid crashes
            c_im_found = .false.
            return
        end if
    end if
    !if we got here it means that we're looking at a corner solution (which the current version of the
    !program does not handle (and the program only gets to this point if we supply wrong parameters).

    write(*,*) 'Error in subroutine get_c_im : not generalized for corner solutions in taxhetgrtr.'
    sync all
    error stop



end subroutine get_c_im

!subroutine for generating grids. 4 equispaced grids are generated for the first four states
!(a_1,a_2,rho,theta).
!The grid for shock space is loaded from parameters and we keep it separate (this grid is
!never resized).
subroutine grids_gen(A1_gr,rho_gr,T_gr,N_a,N_rho,N_T,pars)
    real(dp), dimension(:), intent(out) :: A1_gr,rho_gr,T_gr
    integer, intent(in) :: N_a,N_rho,N_T !can't just load these from parameter type
    !because in CTS algorithm, they vary.
    type(par) :: pars

    call gen_equi_grid(A1_gr,pars%a_min(1,1),pars%a_max(1,1))
    call gen_equi_grid(rho_gr,pars%rho_min,pars%rho_max)

    !Grid for 'time' is potentially a bit more complicated
    if(N_T==1) then !this should happen only if static_environment == 1.
        T_gr(1) = 0.0_dp !Fix time at t=0.
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
!as to whether we are on the grid (to increase efficiency). Later on can introduce
!some debugging option but should do this using pre-processing (conditional
!compilation).
subroutine interp_V(V,x,s_ind,y,a1_gr,rho_gr,t_gr,N_a,N_rho,N_t)
    real(dp), dimension(:,:,:,:), intent(in) :: V !4-dimensional value function if I=2
    real(dp), dimension(3), intent(in) :: x !3 variables which we interpolate over (no interp over M)
    integer, intent(in) :: s_ind
    real(dp), intent(out) :: y
    real(dp), dimension(:) ,intent(in) :: a1_gr,rho_gr,t_gr

    integer, intent(in) :: N_a,N_rho,N_t
    !Same number of gridpoints for agent 1 and agent 2 is assumed - this can be
    !generalized easily
    integer :: j_a1,j_rho,j_t !used for storing which points are to be used in interpolation.
    real(dp) :: t_a1,t_rho,t_t ! weights in linear combinations t*x(j) + (1-t)*x(j+1)


    !This follows quite closely algorithm suggested in the Numerical recipes in Fortran book.
    !1) Find the points to be used in interpolation (for all variables, point with index j
    !is the point s.t. x lies between x(j) and x(j+1).
    call locate_gr(a1_gr,N_a,x(1),j_a1) !j_a1 will be the index of the grid point in a1_gr
    call locate_gr(rho_gr,N_rho,x(2),j_rho)


    !coefficients in linear combinations (value 0 is associated with the lower point)
    t_a1 = (a1_gr(j_a1+1)-x(1))/(a1_gr(j_a1+1)-a1_gr(j_a1))
    t_rho = (rho_gr(j_rho+1)-x(2))/(rho_gr(j_rho+1)-rho_gr(j_rho))

    !split The algorithm into 2 branches depending on whether we are also
    !interpolating over theta or not!
    if(N_t>1) then
        !3-dimensional interpolation (when I=2, needs to be generalized for I=3)
        call locate_gr(t_gr,N_t,x(3),j_t)
        t_t = (t_gr(j_t+1)-x(3))/(t_gr(j_t+1)-t_gr(j_t))

        y = V(j_a1,j_rho,j_t,s_ind)*t_a1*t_rho*t_t+&
            V(j_a1+1,j_rho,j_t,s_ind)*(1-t_a1)*t_rho*t_t+&
            V(j_a1,j_rho+1,j_t,s_ind)*t_a1*(1-t_rho)*t_t+&
            V(j_a1+1,j_rho+1,j_t,s_ind)*(1-t_a1)*(1-t_rho)*t_t+&
            V(j_a1,j_rho,j_t+1,s_ind)*t_a1*t_rho*(1-t_t)+&
            V(j_a1+1,j_rho,j_t+1,s_ind)*(1-t_a1)*t_rho*(1-t_t)+&
            V(j_a1,j_rho+1,j_t+1,s_ind)*t_a1*(1-t_rho)*(1-t_t)+&
            V(j_a1+1,j_rho+1,j_t+1,s_ind)*(1-t_a1)*(1-t_rho)*(1-t_t)

    else
        !2-dimensional interpolation. In this case, the third index
        !is fixed, so we work with V(:,:,1,s_ind).

        y = V(j_a1,j_rho,1,s_ind)*t_a1*t_rho+&
            V(j_a1+1,j_rho,1,s_ind)*(1-t_a1)*t_rho+&
            V(j_a1,j_rho+1,1,s_ind)*t_a1*(1-t_rho)+&
            V(j_a1+1,j_rho+1,1,s_ind)*(1-t_a1)*(1-t_rho)
    end if

end subroutine interp_V
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
subroutine LFFC(a_st,rho_st,t_st,g_st,g_st_ind,pars,b_min_mat,b_max_mat,l_max_mat,a_min_mat,a_max_mat,&
rho_min_vec,rho_max_vec,LFFC_choices,num_fc,WB_poss,c_guess_all_gp,fc_found,bcsf_loss)
    !Values of state variables as saved in the grid (hence the suffix st).
    real(dp), dimension(I_par-1), intent(in) :: a_st
    real(dp), dimension(I_par-1), intent(in) :: rho_st
    real(dp), intent(in) :: t_st,g_st
    integer, intent(in) :: g_st_ind
    type(par), intent(in) :: pars
    real(dp), dimension(M_par,I_par), intent(in) :: l_max_mat,b_min_mat,b_max_mat
    real(dp), dimension(M_par,I_par-1), intent(in) :: a_min_mat,a_max_mat
    real(dp), dimension(M_par,1), intent(in) :: rho_min_vec,rho_max_vec
    real(dp), dimension(:,:,:), intent(in) :: LFFC_choices
    integer, intent(in) :: num_fc
    integer, dimension(1+I_par*M_par,3) :: WB_poss !possibilities of what_binds
    type(c_guess_all), intent(out) :: c_guess_all_gp !output: all feasible choices at a grid point
    logical, intent(out) :: fc_found !true if a feasible choice found at this gridpoint


    real(dp), intent(out) :: bcsf_loss !minimum loss (that will either be zero if strictly feasible
    !point found or the loss at the minimum loss point that was found).

    real(dp), dimension(1,I_par) :: theta_lvl
    real(dp), dimension(M_par,I_par) :: max_cons_mat
    real(dp), dimension(M_par,I_par) :: g_lvl_mat
    integer :: m_ind
    real(dp), dimension(M_par,M_par) :: P_onerow
    logical :: c_im_found

    integer :: WB_ind,wb_ind_max, fc_ind !what_binds index and index for feasible choice (used to loop through LFFC_choices)

    logical :: constr_ok
    real(dp) :: loss

    !variables used for finding best choice (the least infeasible in the sense of having
    !the lowest loss) (bcsf stands for best choice so far)
    integer :: bcsf_wb_ind !what_binds index
    integer, dimension(3) :: bcsf_wb !what_binds for the best choice so far
    real(dp), dimension(M_par,I_par) :: bcsf_c_gp

    real(dp), dimension(M_par,I_par) :: l,b_prime
    real(dp), dimension(M_par,I_par-1) :: a_prime
    !real(dp), dimension(1,I_par) :: theta_prime
    real(dp), dimension(M_par,I_par-1) :: rho_prime

    !Local choice variables: gp stands for grid point.
    real(dp), dimension(M_par,I_par) :: c_gp



    !Start by computing maximum feasible consumption using max_feas_cons -> this will be the same
    !for all choices at this gridpoint and will be used to compute the actual choices.
    !As prerequisites we need productivity in levels and government expenditure in levels
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

    !If we look for interior solutions only, only check the first WB possibility which corresponds to
    !interior solution
    if(pars%interior_solution_only) then
        wb_ind_max = 1
    else
        wb_ind_max = 1+I_par*M_par
    end if


    do WB_ind = 1,wb_ind_max
        do fc_ind = 1,num_fc !Try all consumption choices in LFFC_choices
            c_gp = LFFC_choices(:,:,fc_ind)*max_cons_mat
            !c_gp contains the actual consumption choices by agents (corresponding to c's in the paper).


            !We need to get the other variables. First of all get the last element of the
            !consumption matrix. This depends on what_binds.
            c_im_found = .false.



            call get_c_im(c_gp,WB_poss(WB_ind,:),a_st,rho_st,theta_lvl,g_lvl_mat(:,1)&
            ,g_st_ind,pars,P_onerow,c_im_found)


            !c_im_found will be .false. if the iterative algorithm couldn't find c from constraint
            !in this case skip to the next index fc_ind
            if(c_im_found == .false.) go to 100

            !Should have one subroutine which, for a given consumption choice,
            !computes all the remaining choice variables and next-period states.
            !The following bits should later be put in such a subroutine and it should be used
            !everywhere in the program to ensure consistency!!!
            !The input should be the complete consumption matrix because we might want to treat
            !exceptions (such as negative c_IM or failure to find root of nonlinear equation)
            !differently depending on the place!!!

            !If we got here we have a full consumption matrix in c_gp. Compute all the
            !remaining variables and check feasibility.
            call lab(c_gp,l,theta_lvl,pars%mass,g_lvl_mat(:,1),pars%Gamma_par)


            call ass(a_st,c_gp,l,pars%A_par,pars%B_par,pars%Gamma_par,pars%beta,P_onerow,a_prime)


            !getting assets in levels at this point will become necessary only once
            !borrowing constraints extension is introduced. Otherwise we can work with a's.

            !rho_prime
            call get_rho_prime(c_gp,rho_prime)


            call check_constr(c_gp,l,l_max_mat,theta_lvl,a_prime,a_min_mat&
            ,a_max_mat,WB_poss(WB_ind,:),rho_prime,rho_min_vec,rho_max_vec,constr_ok,loss,pars%A_par,fc_found)
            !fc_found: as long as it is false, we go through the whole function check_constr (computing
            !loss) because we may need to use a choice where some constraint is slightly violated.

            !If the choice is feasible, save it and go to 200 (try finding feasible choices for next point)
            if(constr_ok) then
                fc_found = .true. !we found at least one feasible choice
                !Save the result into c_guess_all_gp which contains all feasible guesses
                !at a gridpoint (or one as little infeasible guess as possible).
                c_guess_all_gp%c_guess_all(wb_ind)%c_guess = c_gp
                c_guess_all_gp%c_guess_all(wb_ind)%what_binds = WB_poss(WB_ind,:)
                bcsf_loss = loss !this should be zero

                !If we were only looking for interior solution or if we are satisfied with
                !one feasible choice per gridpoint, return
                if(pars%one_fc_suff) return
                !go to next what_binds possibility
                go to 200

            else !constr_ok == .false., so choice not strictly feasible.
            !However, if we didn't find a feasible choice yet we may want to keep this, as
            !we want LFFC to return the least infeasible choice if no feasible choice was found.

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


            end if

100     end do
        !If looking for interior solution only, skip this step -> I.e. break cycle or exit or whatever it's called.
200 end do


!If no feasible choice found, save the least infeasible choice and return loss so we can compute average loss
!and get an idea of how bad the violation of constraints is
if(.not. fc_found) then
c_guess_all_gp%c_guess_all(bcsf_wb_ind)%c_guess = bcsf_c_gp
c_guess_all_gp%c_guess_all(bcsf_wb_ind)%what_binds = bcsf_wb
end if


end subroutine LFFC
!
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


!What we worked with until now are shares of consumption of all (potentially) available resources
!consumed by a group of agents (so that c_i*pi_i = something). It is more convenient to save these
!choices already normalized by mass of agents, so we can get the actual consumption choices
!simply by multiplying the saved share by the maximum consumption, and we won't have to divide by mass
!every time.
!Do this for every agent
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

!subroutine load_V loads the value function (in vector form) and grids from files.
subroutine load_V(V_vect,a1_gr,rho_gr,t_gr,input_folder,folder_name)
    real(dp), dimension(:), intent(out) :: V_vect
    real(dp), dimension(:), intent(out) :: a1_gr,rho_gr,t_gr

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

    file_path = 'results/'//trim(input_folder)//'/rho_gr.out'
    open(unit=20, file = file_path, status = 'old')
    read(unit=20, fmt=*)  rho_gr
    close(unit = 20)

    file_path = 'results/'//trim(input_folder)//'/t_gr.out'
    open(unit=20, file = file_path, status = 'old')
    read(unit=20, fmt=*)  t_gr
    close(unit = 20)

end subroutine load_V


!Subroutine max_feas_cons computes maximum feasible consumption going to a group
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
!
!The first element of a_min and a_max is the
!constraint on a2. Rho is a column vector (needs
!to be generalized to matrix when I>2). Also, borrowing constraints are not yet implemented in the model with transfers
!but the bounds on them are left here (their checking will be left out of checking constraints).
subroutine prep_bounds(l_max,b_min,b_max,a_min,a_max,rho_min,rho_max,l_max_mat,b_min_mat,b_max_mat,&
a_min_mat,a_max_mat,rho_min_vec,rho_max_vec)
    real(dp), dimension(1,I_par), intent(in) :: l_max,b_min,b_max
    real(dp), dimension(1,I_par-1), intent(in) :: a_min,a_max
    real(dp), intent(in) :: rho_min,rho_max

    real(dp), dimension(M_par,I_par), intent(out) :: l_max_mat,b_min_mat,b_max_mat
    real(dp), dimension(M_par,I_par-1), intent(out) :: a_min_mat,a_max_mat
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
    do i=1,I_par-1
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

!Function t_prime is a trivial function capturing the law of motion of 'time' state. t is incremented by one
!unless it would exceed t_max, in which case it will be set equal to t_max.
function t_prime(t,t_max)
    real(dp) :: t_prime
    real(dp), intent(in) :: t,t_max

    t_prime = minval([t + 1.0_dp,t_max])

end function t_prime


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
    real(dp), dimension(:,:,:,:), intent(in) :: V_old !4 states if I=2
    type(grids_type), intent(in) :: grids_old,grids_new

    real(dp), dimension(:,:,:,:), intent(out) :: V_new

    character*(*), intent(in) :: folder_name

    integer, dimension(4) :: shape_old,shape_new,shape_grid_old,shape_grid_new

    integer :: a1_ind,a2_ind,rho_ind,t_ind,g_ind

!Perform checks that the grids and the value functions are compatible in the sense of dimensionality.


    shape_old = shape(V_old)
    shape_new = shape(V_new)

    shape_grid_old(1) = size(grids_old%a1_gr)
    shape_grid_old(2) = size(grids_old%rho_gr)
    shape_grid_old(3) = size(grids_old%t_gr)
    !(grids for shock realization not checked). They should always
    !be the same and we only check the whether the value functions
    !have the same M_par. Also, no interpolation happens over shock
    !realizations - we just assume that shock space is the same here.

    shape_grid_new(1) = size(grids_new%a1_gr)
    shape_grid_new(2) = size(grids_new%rho_gr)
    shape_grid_new(3) = size(grids_new%t_gr)

    if(any(shape_old(1:3)/=shape_grid_old(1:3))) then
        call double_output(folder_name, 'Error: grids_old in subroutine V_from_V incompatible with V_old.')
        error stop
    end if
    if(any(shape_new(1:3)/=shape_grid_new(1:3))) then
        call double_output(folder_name, 'Error: grids_new in subroutine V_from_V incompatible with V_new.')
        error stop
    end if

    !If any of the gridpoints has different M_par, give an error. These should never happen but check
    !just in case, because a bug here could be difficult to detect.
    if(shape_old(4) /= M_par) then
        call double_output(folder_name, 'Error: V_old in V_from_V computed for different M_par than the current one.')
        error stop
    end if
    !If any of the gridpoints has different M_par, give an error.
    if(shape_new(4) /= M_par) then
        call double_output(folder_name, 'Error: V_new in V_from_V computed for different M_par than the current one.')
        error stop
    end if
    !ALSO check whether the ranges of the new grid exceed those of the old grid. If yes,
    !extrapolation is used. For small differences this is fine but a warning should be displayed.
    !(do not assume monotonicy in grids)
    if((minval(grids_new%a1_gr)<minval(grids_old%a1_gr)) .or. &
    (minval(grids_new%rho_gr)<minval(grids_old%rho_gr)) .or. &
    (minval(grids_new%t_gr)<minval(grids_old%t_gr)) .or. &
    (maxval(grids_new%a1_gr)>maxval(grids_old%a1_gr)) .or. &
    (maxval(grids_new%rho_gr)>maxval(grids_old%rho_gr)) .or. &
    (maxval(grids_new%t_gr)>maxval(grids_old%t_gr))) then
        call double_output(folder_name,'Warning: in V_from_V, bounds of grid_new exceed bounds of grid_old.')
        call double_output(folder_name,'Extrapolation will be used, dangerous far from bounds.')
    end if



    do a1_ind = 1,shape_grid_new(1)
        do rho_ind = 1,shape_grid_new(2)
            do t_ind = 1,shape_grid_new(3)
                do g_ind = 1,M_par
                    !get the interpolated value and save it in the appropriate element of V_new
                    call interp_V(V_old,[grids_new%a1_gr(a1_ind),grids_new%rho_gr(rho_ind),&
                        grids_new%t_gr(t_ind)],g_ind,V_new(a1_ind,rho_ind,t_ind,g_ind)&
                        ,grids_old%a1_gr,grids_old%rho_gr,grids_old%t_gr,&
                        shape_grid_old(1),shape_grid_old(2),shape_grid_old(3))
                        !interp_V(V,x,s_ind,y,a1_gr,rho_gr,t_gr,N_a,N_rho,N_t)
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

    integer, dimension(4), intent(in) :: ind_unfold !vector of indices (this assumes I=2, generalization will be needed)
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
    real(dp), dimension(:,:,:,:), pointer :: V_old_pntr

    !optimal index of what_binds (this is used in Howard acceleration algorithm)
    integer, intent(inout) :: optim_wb_gp

    !Values of state variables as saved in the grid (hence the suffix st).
    real(dp), dimension(I_par-1) :: a_st
    real(dp), dimension(I_par-1) :: rho_st
    real(dp) :: t_st,g_st
    integer :: g_st_ind
    !prepared bounds...
    real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat
    real(dp), dimension(M_par,I_par-1) :: a_min_mat,a_max_mat
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

    !real(dp), dimension(M_par,I_par) :: c_tmp_debug

    !Common block to be accessible by subroutine CPC_ret
    common /cpc_ret_bloc/ a_st,theta_lvl,rho_st,t_st,g_st,g_st_ind
    !(EFE) : The bounds could actually be prepared even before this subroutine is called,
    !they are the same at every gridpoint.
    common /cpc_ret_bloc2/ g_lvl_mat, P_onerow,l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec

    common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
    !to some subroutines.

    common /commloss/ loss


    !If the number of free variables in the optimization problem is 1, then
    !we can't use subroutine e04jcf. This will be the case only when we are solving
    !deterministic version of the model with two players.
    if(M_par*I_par-1 == 1) then
        optim_subroutine = 3
    end if

    !Much of what follows is similar to what was done elsewhere (LFFC).

    !Get current state from grid:
    a_st(1) = grids%a1_gr(ind_unfold(1))
    rho_st = grids%rho_gr(ind_unfold(2))
    t_st = grids%t_gr(ind_unfold(3))
    g_st = pars%S(1,ind_unfold(4))
    g_st_ind = ind_unfold(4)


    !prepare bounds for checking constraints
    call prep_bounds(pars%l_max,pars%b_min,pars%b_max,pars%a_min,pars%a_max,pars%rho_min,pars%rho_max,&
    l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec)
    !constraints on b are not even used at this point but they might be in the future, therefore
    !b_m**_mat in from previous versions.


    call prod_lvl(t_st,pars%theta_0,pars%xi_theta,pars%static_environment,theta_lvl)

    !Need government expenditure in levels for all shock realizations. Actually save it in a MxI matrix
    !where every row contains gov. exp. in level in state m
    do m_ind = 1,M_par
        g_lvl_mat(m_ind,:) = pars%S(1,m_ind)*trend_g(pars%mass,theta_lvl,pars%l_max,pars%k_gov)
    end do


!    !P_onerow = matmul(reshape([1.0_dp,1.0_dp],[M_par,1]),reshape(pars%P(g_st_ind,:),[1,M_par]))
    P_onerow = pars%P_onerow(g_st_ind,:,:)



    !If we are skipping maximization at this point, compute the 'optimal' value using previous choice
    !and return (Howard).
    if(skip_max) then

        write(*,*) 'Error: skip_max not implemented yet in VFI_solve_GP!!!'
        !When implementing it, look at TaxHetGr program
        error stop

    end if


    call max_feas_cons(pars%l_max,pars%mass,theta_lvl,g_lvl_mat,max_cons_mat)
    !Use this to generate the vector of upper bounds...

    !get upper bounds
    call c_inc_mtv(max_cons_mat,bu)

    V_new_gp = -10.0_dp**20 !initialize value (must be very low to make sure that this
    !does not serve as an artificial lower bound and also so we detect errors
    !arising from its useage easily). Also if at some point we didn't
    !find a value better than this, then the optimal policy would not be updated.


    if(glob_opt) then !use global optimization
        !Call initialization routing (has to be done only once)
        call e05jaf(0,comm,lcomm,ifail)
    end if

    !Need to cycle over elements of what_binds. Get a maximum at every index.
    !If there is -1 (no guess found for that what_binds), skip to the next row of what_binds...
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

              write(*,*) 'Error in VFI_solve_gp: CPC_ret3 not yet implemented'
              error stop
              !call e05jbf(M_par*I_par-1,CPC_ret3,0,0,bl,bu,sdlist,list,numpts,initpt,E05JBK, &
              !  xc,fmax,comm,lcomm,iuser,ruser,ifail)

                !if ifail>0, then fmax could contain some nonsense (and xc could be nonsense).
                !If this happens, perhaps use one of the local minimization subroutines
                !starting at the initial guess...
        else

        select case(optim_subroutine)
            case(1) !Quadratic approximation method
!              write(*,*) 'Error in VFI_solve_gp: CPC_ret2 not yet implemented'
!              error stop
                call e04jcf(CPC_ret2,M_par*I_par-1,ncvar_par*2+1,xc,bl,bu,0.1_dp,&
                0.000001_dp,e04jcp,200,fmax,nf,iuser,ruser,ifail)
            case(2) !Newton-like method
                call e04jyf(M_par*I_par-1,0,CPC_ret,bl,bu,xc,fmax,iw,liw,w,lw,iuser,ruser,ifail)
            case(3) !1-D optimization (only used then M_par = 1, I_par = 2)
                A = bl(1)
                B = bu(1)
                e1 = 0.000000001
                e2 = 0.0000000001
                maxcal = 100

                call E04ABA(CPC_ret_1d,E1,E2,A,B,MAXCAL,XC(1),fmax,IUSER,RUSER,IFAIL)

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

                !Save optimal wb index (optimal in the sense that even though maximization
                !failed to find something good, initial guess for this grid point was the
                !best).
                optim_wb_gp = guess_ind


            end if

        end if

100 end do



end subroutine VFI_solve_gp



!Subroutine solve_FP_problem solves the first period problem (t=0), given a value function
!for periods t>=1 and initial conditions.
!
!It is better to have one subroutine for solving the FP problem and one for solving
!the subsequent periods because the dimensionality of choice variables is different
!and there are a lot of other small changes which would make it confusing if we
!put all simulation into one subroutine).
subroutine solve_FP_problem(V_old,b_init,theta_0,pars,FP_result,t_gr,g_ind_init)
    Use nag_library, Only: e04jyf, nag_wp, e04jcf, e04jcp,e05jaf,e05jbf,E05JBK

    real(dp), dimension(:,:,:,:), intent(in) :: V_old
    real(dp), dimension(1,I_par), intent(in) :: b_init
    real(dp), dimension(1,I_par), intent(in) :: theta_0 !initial productivity (in levels)
    type(par), intent(in) :: pars
    type(FP_result_type), intent(out) :: FP_result
    real(dp), dimension(:), intent(in) :: t_gr

    integer, intent(in) :: g_ind_init

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
    real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat
    real(dp), dimension(M_par,I_par-1) :: a_min_mat,a_max_mat
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



    !Debug - possibly remove later
    real(dp), dimension(1,I_par) :: c_ss,l_ss

    !Debug values
    real(dp) :: fc_debug
    real(dp), dimension(I_par) :: xc_debug
    integer :: inform_debug

    !Net asset positions
    real(dp), dimension(1,I_par-1) :: b_net_init,b_net_prime


    !Declare the common block use in CPC return computation.
    common /cpc_ret_fp/ b_init_cp, theta_0_cp,g_ind_init_cp

    !common /cpc_ret_fp_bnd/ l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec
    common /cpc_ret_bloc_fp/ g_lvl_mat,l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec

    !Create copies of inputs which need to be made available to subroutine
    !called by the maximization subroutine using a common block.
    b_init_cp = b_init
    theta_0_cp = theta_0
    g_ind_init_cp = g_ind_init

    !If g_ind_init = -1, then draw g_ind_init randomly corresponding to probability in a stationary distribution.
    !Otherwise just set it to input g_ind_init



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


    !Initialization routine for e05jbf
    call e05jaf(0,comm,lcomm,ifail)

    !Need to have bounds prepared...... And call the initialization subroutine...
    ifail = -1 !set 1 for silent exit
    numpts = 100



    !Run global optimization subroutine with lots of different initial conditions
    !(they are random) - we can afford this robustness in the first period
    fmax_tmp = 100000000000000 !initialize (working with minus 1 multiples so lower number is better here)
    xc_tmp = 0.0000000001_dp
    do max_ind = 1,pars%sim_attmax1 !This is possibly quite time-consuming.
        ifail = 1
        call e05jbf(I_par,CPC_ret3_FP,0,4,bl,bu,sdlist,list,numpts,initpt,E05JBK, &
            xc,fmax,comm,lcomm,iuser,ruser,ifail)
        if(fmax<fmax_tmp) then
            fmax_tmp = fmax
            xc_tmp = xc

            !debug
            !write(*,*) 'fmax=',fmax,'c=', xc
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

!write(*,*) 'fmax2=',fmax2,'c=', xc2

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


    !get the initial net asset position (t=-1) from the gross ones:
    b_net_init(1,1) = b_init(1,2) - b_init(1,1)

    !Get next-period asset holdings:
    call b_prime_fp(b_net_init,FP_result%c,FP_result%l,pars%A_par,pars%B_par,pars%Gamma_par,FP_result%b_net_prime)

    !Get marginal utility of consumption
    call util_c(FP_result%c,FP_result%u_c,pars%A_par)
    !get a_prime (for I=2, needs to be generalized for I=3)
    FP_result%a_prime = (FP_result%b_net_prime*FP_result%u_c(1,2))/pars%beta

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
        integer, intent(in) :: nstate !=1 if the subroutine is called for the first time,
        integer, intent (in) :: n
        real (dp), intent (inout) :: ruser(*)
        real (dp), intent (in) :: xc(n)
        integer, intent (inout) :: iuser(*) !In this instance, it is used to pass what_binds
        integer, intent(out) :: inform


        real(dp), dimension(:,:,:,:), pointer :: V_old_pntr !pointer to V_old, passed using common block

        !The following variables are passed to the function using a common block
        real(dp), dimension(:), pointer :: a1_gr_pntr,rho_gr_pntr,t_gr_pntr
        real(dp), dimension(:,:), pointer :: S_pntr,P_pntr, mass_pntr
        !and transition matrix.
        real(dp), pointer :: A_par_pntr,B_par_pntr,Gamma_par_pntr
        type(par), pointer :: pars_pntr

        !Values of state variables as saved in the grid (hence the suffix st).
        real(dp), dimension(I_par-1) :: a_st
        real(dp), dimension(I_par-1) :: rho_st
        real(dp) :: t_st,g_st
        integer :: g_st_ind
        real(dp), dimension(1,I_par) :: theta_lvl

        real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat
        real(dp), dimension(M_par,I_par-1) :: a_min_mat,a_max_mat
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
        real(dp), dimension(1,I_par) :: b_prime
        real(dp), dimension(1,I_par-1) :: a_prime
        real(dp), dimension(M_par,1) :: rho_prime !with more states it's a column vector

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

        !Net asset positions
        real(dp), dimension(1,I_par-1) :: b_net_init,b_net_prime



        common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
        !to some subroutines.
        common /pars_pointer/ pars_pntr

        common /cpc_ret_bloc/ a_st,theta_lvl,rho_st,t_st,g_st,g_st_ind
        common /cpc_ret_bloc_fp/ g_lvl_mat,l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec


        common /comm1/ a1_gr_pntr,rho_gr_pntr,t_gr_pntr,A_par_pntr,B_par_pntr,&
        Gamma_par_pntr,S_pntr,P_pntr,mass_pntr

        common /cpc_ret_fp/ b_init_cp, theta_0_cp,g_ind_init_cp

        !Negative value lead to termination of maximization procedure.
        inform = 0

        c_vec = reshape(xc,[1,I_par]) !Convenient to have this so we don't have to use
        !reshape all the time (element (1,i) contains consumption of agent i in period 0)

        !(Then generalizing law of motion for productivity, the producivity in first period
        !will be for index t=0).

        !Given consumption, compute labour supply (using resource constraint and intratemporal EE)
        call lab_FP(c_vec,lab_vec,theta_0_cp,mass_pntr,g_lvl_mat(g_ind_init_cp,1),Gamma_par_pntr)

        !get the initial net asset position (t=-1) from the gross ones:
        b_net_init(1,1) = b_init_cp(1,2) - b_init_cp(1,1)

        !Get next-period asset holdings (now b_prime are already the net asset holdings!)
        call b_prime_fp(b_net_init,c_vec,lab_vec,A_par_pntr,B_par_pntr,Gamma_par_pntr,b_net_prime)


        !Get marginal utility of consumption
        call util_c(c_vec,u_c_vec,A_par_pntr)


        if(I_par>2) then
            write(*,*) 'Error: subroutines CPC_ret_fp** need to be generalized for I>2 case.'
            error stop
        end if

        !get a_prime
        a_prime = (b_net_prime*u_c_vec(1,2))/pars_pntr%beta

        !This will need to be generalized for I>2
        rho_prime(1,1) = u_c_vec(1,1)/u_c_vec(1,2)


        !Get next-period 'time' (assumes that the first point of the grid is the time in the
        !second period (index t=1)
        t_st_pr = t_gr_pntr(1)

        !utility
        call util(c_vec,lab_vec,U,A_par_pntr,B_par_pntr,Gamma_par_pntr) !need to use mod_util for this


        !Current return
        fc = sum(pars_pntr%mass*pars_pntr%alpha*U)




        call check_constr(matmul(ones,c_vec),matmul(ones,lab_vec)&
        ,l_max_mat,theta_0_cp,matmul(ones,a_prime),a_min_mat&
        ,a_max_mat,iuser,rho_prime,rho_min_vec,rho_max_vec,constr_ok,loss,pars_pntr%A_par,.false.)




        !Divide loss by M_par (so it's not counted here more than in later periods)
        fc = fc - loss/real(M_par,dp)


        !2) Truncate (instead of extrapolation, use nearest neighbour extrap + loss)
        !(just like in CPC_ret and other functions but have to use reshape here
        !because the constraints are in matrix form

        if(.not. constr_ok) then
            a_prime=max(reshape(a_min_mat(1,:),[1,1]),a_prime)
            a_prime=min(reshape(a_max_mat(1,:),[1,1]),a_prime)
            rho_prime=max(reshape(rho_min_vec(1,:),[M_par,1]),rho_prime)
            rho_prime=min(reshape(rho_max_vec(1,:),[M_par,1]),rho_prime)
        end if


        !3) Use interpolation to get the continuation value
        call interp_V(V_old_pntr,[a_prime(1,1),rho_prime(1,1),t_st_pr]&
        ,g_ind_init_cp,conval,a1_gr_pntr,rho_gr_pntr,t_gr_pntr,size(a1_gr_pntr),&
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


        real(dp), dimension(:,:,:,:), pointer :: V_old_pntr !pointer to V_old, passed using common block

        !The following variables are passed to the function using a common block
        real(dp), dimension(:), pointer :: a1_gr_pntr,rho_gr_pntr,t_gr_pntr
        real(dp), dimension(:,:), pointer :: S_pntr,P_pntr, mass_pntr
        !and transition matrix.
        real(dp), pointer :: A_par_pntr,B_par_pntr,Gamma_par_pntr
        type(par), pointer :: pars_pntr

        !Values of state variables as saved in the grid (hence the suffix st).
        real(dp), dimension(I_par-1) :: a_st
        real(dp), dimension(I_par-1) :: rho_st
        real(dp) :: t_st,g_st
        integer :: g_st_ind
        real(dp), dimension(1,I_par) :: theta_lvl

        real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat
        real(dp), dimension(M_par,I_par-1) :: a_min_mat,a_max_mat
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
        real(dp), dimension(1,I_par) :: b_prime
        real(dp), dimension(1,I_par-1) :: a_prime
        real(dp), dimension(M_par,1) :: rho_prime !with more states it's a column vector

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

        !Net asset positions
        real(dp), dimension(1,I_par-1) :: b_net_init,b_net_prime



        common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
        !to some subroutines.
        common /pars_pointer/ pars_pntr

        common /cpc_ret_bloc/ a_st,theta_lvl,rho_st,t_st,g_st,g_st_ind
        common /cpc_ret_bloc_fp/ g_lvl_mat,l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec


        common /comm1/ a1_gr_pntr,rho_gr_pntr,t_gr_pntr,A_par_pntr,B_par_pntr,&
        Gamma_par_pntr,S_pntr,P_pntr,mass_pntr

        common /cpc_ret_fp/ b_init_cp, theta_0_cp,g_ind_init_cp

        !Negative value lead to termination of maximization procedure.
        inform = 0

        c_vec = reshape(xc,[1,I_par]) !Convenient to have this so we don't have to use
        !reshape all the time (element (1,i) contains consumption of agent i in period 0)

        !(Then generalizing law of motion for productivity, the producivity in first period
        !will be for index t=0).

        !Given consumption, compute labour supply (using resource constraint and intratemporal EE)
        call lab_FP(c_vec,lab_vec,theta_0_cp,mass_pntr,g_lvl_mat(g_ind_init_cp,1),Gamma_par_pntr)

        !get the initial net asset position (t=-1) from the gross ones:
        b_net_init(1,1) = b_init_cp(1,2) - b_init_cp(1,1)

        !Get next-period asset holdings (now b_prime are already the net asset holdings!)
        call b_prime_fp(b_net_init,c_vec,lab_vec,A_par_pntr,B_par_pntr,Gamma_par_pntr,b_net_prime)


        !Get marginal utility of consumption
        call util_c(c_vec,u_c_vec,A_par_pntr)


        if(I_par>2) then
            write(*,*) 'Error: subroutines CPC_ret_fp** need to be generalized for I>2 case.'
            error stop
        end if

        !get a_prime
        a_prime = (b_net_prime*u_c_vec(1,2))/pars_pntr%beta

        !This will need to be generalized for I>2
        rho_prime(1,1) = u_c_vec(1,1)/u_c_vec(1,2)


        !Get next-period 'time' (assumes that the first point of the grid is the time in the
        !second period (index t=1)
        t_st_pr = t_gr_pntr(1)

        !utility
        call util(c_vec,lab_vec,U,A_par_pntr,B_par_pntr,Gamma_par_pntr) !need to use mod_util for this


        !Current return
        fc = sum(pars_pntr%mass*pars_pntr%alpha*U)


        !Continuation return
        !1) Get the loss associated with violating constraints



        call check_constr(matmul(ones,c_vec),matmul(ones,lab_vec)&
        ,l_max_mat,theta_0_cp,matmul(ones,a_prime),a_min_mat&
        ,a_max_mat,iuser,rho_prime,rho_min_vec,rho_max_vec,constr_ok,loss,pars_pntr%A_par,.false.)

        
        fc = fc - loss/real(M_par,dp)


        !2) Truncate (instead of extrapolation, use nearest neighbour extrap + loss)
        !(just like in CPC_ret and other functions but have to use reshape here
        !because the constraints are in matrix form

        if(.not. constr_ok) then
            a_prime=max(reshape(a_min_mat(1,:),[1,1]),a_prime)
            a_prime=min(reshape(a_max_mat(1,:),[1,1]),a_prime)
            rho_prime=max(reshape(rho_min_vec(1,:),[M_par,1]),rho_prime)
            rho_prime=min(reshape(rho_max_vec(1,:),[M_par,1]),rho_prime)
        end if



        !3) Use interpolation to get the continuation value
        call interp_V(V_old_pntr,[a_prime(1,1),rho_prime(1,1),t_st_pr]&
        ,g_ind_init_cp,conval,a1_gr_pntr,rho_gr_pntr,t_gr_pntr,size(a1_gr_pntr),&
        size(rho_gr_pntr),size(t_gr_pntr))



        fc = fc + pars_pntr%beta*conval!+ discounted continuation value (loss vas already subtracted and is not discounted)


        !At the end multiply fc by -1 (because we want to use it in a minimization subroutine)
        fc = fc*(-1.0_dp)


end subroutine CPC_ret2_FP

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
    !that will be passed between periods.

    real(dp), dimension(M_par,I_par) :: c_guess_old !old c_guess (used in cases where maximization
    !yields nonsensical result (such as negative consumption in last state of agent I)

    real(dp), dimension(M_par,1) :: ones_col = 1.0_dp

    integer :: t_ind !period index

    !Pointer to value function. This will be used to make the old value function V_old
    !accessible to the subroutine used in NAG maximization.
    real(dp), dimension(:,:,:,:), pointer :: V_old_pntr


    !Values of state variables as saved in the grid (hence the suffix st).
    real(dp), dimension(I_par-1) :: a_st
    real(dp), dimension(I_par-1) :: rho_st
    real(dp) :: t_st,g_st,t_st_pr
    integer :: g_st_ind
    !prepared bounds...
    real(dp), dimension(M_par,I_par) :: l_max_mat,b_min_mat,b_max_mat
    real(dp), dimension(M_par,I_par-1) :: a_min_mat,a_max_mat
    real(dp), dimension(M_par,1) :: rho_min_vec,rho_max_vec

    real(dp), dimension(1,I_par) :: theta_lvl
    real(dp), dimension(M_par,I_par) ::  g_lvl_mat,max_cons_mat
    real(dp) :: theta_st_pr !next-period relative productivitiy of agent 2
    real(dp), dimension(M_par,M_par) :: P_onerow

    !guess_ind used to cycle over elements of guess_all_gp
    !integer :: guess_ind,best_guess_ind

    integer :: m_ind !for cycling over states

    !Variables used in NAG subroutine E04JYF
    real(dp), dimension(ncvar_par) :: xc,xc_tmp !choice variables (M_par x I_par matrix in row form)
    real(dp) :: fc !function value in maximization
    !The following two arrays are arrays of reals and integers which can be used to supply information
    !to the maximized function as an alternative to COMMON block.
    real(dp), dimension(1) :: ruser !Contains no information, only defined for purposes of using a NAG subroutine
    integer, dimension(3) :: iuser !will be used to pass what_binds to some subroutines.
    integer :: ifail
    !Bounds:
    real(dp), dimension(M_par*I_par-1) :: bl = 0.0001_dp,bu
    real(dp) :: fmax,fmax_tmp
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

    !Current_period choices
    real(dp), dimension(M_par,I_par) :: lab_mat,u_c,u_l
    real(dp), dimension(M_par,I_par-1) :: a_prime,b_net_prime
    real(dp), dimension(M_par,1) :: rho_prime,R

    real(dp), dimension(1,I_par) :: y_agent !product of agent
    real(dp), dimension(1,I_par) :: allinc_agent !all income of agent (labour + net income from assets)

    integer :: max_ind

    !Common block to be accessible by subroutine CPC_ret
    common /cpc_ret_bloc/ a_st,theta_lvl,rho_st,t_st,g_st,g_st_ind

    !(EFE) : The bounds could actually be prepared even before this subroutine is called,
    !they are the same at every gridpoint.
    common /cpc_ret_bloc2/ g_lvl_mat, P_onerow,l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec

    common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
    !to some subroutines.

    !If the number of free variables in the optimization problem is 1, then
    !we can't use subroutine e04jcf. This will be the case only when we are solving
    !deterministic version of the model. with two players

    if(M_par*I_par-1 == 1) then
        optim_subroutine = 2
    end if

!_____________________________________________________________________________
    !Save the first-period results computed previously as first elements of the simulated series
    sim_result%c_sim(1,:) = FP_result%c(1,:)
    sim_result%u_c_sim(1,:) = FP_result%u_c(1,:)
    sim_result%u_l_sim(1,:) = FP_result%u_l(1,:)
    sim_result%l_sim(1,:) = FP_result%l(1,:)
    sim_result%b_net_prime_sim(1,:) = FP_result%b_net_prime(1,:)
    sim_result%a_prime_sim(1,:) = FP_result%a_prime(1,:)
    sim_result%rho_prime_sim(1) = FP_result%rho_prime
    sim_result%t_prime_sim(1) = FP_result%t_prime


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


    !prepare bounds for checking constraints (bounds on b are redundant at this point)
    call prep_bounds(pars%l_max,pars%b_min,pars%b_max,pars%a_min,pars%a_max,pars%rho_min,pars%rho_max,&
    l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec)


    !In periods (indexed from t=2) to T_sim_max, use c_guess from last-period.
    do t_ind=2,pars%T_sim
        !back up guess of consumption (this will either be FP solution or solution from
        !the last period
        c_guess_old = c_guess

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


        P_onerow = pars%P_onerow(g_st_ind,:,:)

        !Transform the initial guess into vector form so that it can be used
        !by NAG subroutines.
        call c_inc_mtv(c_guess,xc)

        iuser = 0 !nothing binds, interior solution. If we want to generalize to corner solutions, then
        !we would need to cycle over what_binds as in VFI_solve_gp. (this variable
        !is used to pass what_binds to some subroutines called by optimizaiton subroutines.


        !run the optimization a couple of times, each time using as the initial guess the result found previously...
        !We can do this here (and don't do it in VFI) because the length of simulated series is not too high
        fmax_tmp = 100000000000000.0_dp !initialize (working with minus 1 multiples so lower number is better here)
        xc_tmp = 0.0000000001_dp
        do max_ind = 1,pars%sim_attmax2
            ifail = 1 !silent exit would be 1
            select case(optim_subroutine)
                case(1) !Quadratic approximation method
                                    call e04jcf(CPC_ret2,M_par*I_par-1,ncvar_par*2+1,xc,bl,bu,0.1_dp,&
                                        0.000001_dp,e04jcp,pars%sim_maxiter,fmax,nf,iuser,ruser,ifail)
                    !write(*,*) 'Error: optim_subroutine = 1 not implemented in sim_series!'
                    !error stop
                case(2) !Newton-like method
                    call e04jyf(M_par*I_par-1,0,CPC_ret,bl,bu,xc,fmax,iw,liw,w,lw,iuser,ruser,ifail)
                case default
                            write(*,*) 'Error in subroutine sim_series - wrong value of optim_subroutine'
                            error stop
            end select
            if(fmax<fmax_tmp) then
                fmax_tmp = fmax
                xc_tmp = xc
                !debug
                !write(*,*) 'fmax=',fmax,'c=', xc
            end if
        end do
        xc = xc_tmp
        fmax = fmax_tmp
        !now xc and fmax is the best choice and function value found by successive minimizations


        c_im_found = .false.
        !Now we have the solution in xc. Need to get the full consumption plan and save
        !it in c_guess so it is used the next period. Assume interior solution here.
        !First get the incomplete consumption matrix (element M_par,I_par missing)
        call c_inc_vtm(xc,c_guess)
        !(iuser contains what_binds, in this case one for interior solution)
        iuser = 0 !just in case some NAg subroutine overwrites it


        call get_c_im(c_guess,iuser,&
            a_st,rho_st,theta_lvl,g_lvl_mat(:,1),g_st_ind,pars,P_onerow,c_im_found)


        !If nothing found, we have a serious problem
        if(.not. c_im_found) then
            write(*,*) 'Warning: get_c_im found no solutions so we cannot get a full consumption plan to be saved.'
            write(*,*) 'Last-period consumption plan used instead.'
            c_guess = c_guess_old
            !c_guess = 0.00001_dp
            !It might be a better idea to reoptimize in this case? Add a loss in return
            !for negative consumption that would be continuous?
            !The initial guess might also yield negative c_{I,M}

            write(*,*) 'Check subroutine sim_series.'
            !error stop
        end if
        !Having the consumption plan (in c_guess), we need to recover all variables of interest.


        call lab(c_guess,lab_mat,theta_lvl,pars%mass,g_lvl_mat(:,1),pars%Gamma_par)


        !Compute next-period states
        t_st_pr = t_prime(t_st,grids%t_gr(pars%N_t))

        call ass(a_st,c_guess,lab_mat,pars%A_par,pars%B_par,pars%Gamma_par,&
        pars%beta,P_onerow,a_prime)

        call get_rho_prime(c_guess,rho_prime)


        call get_b_prime(b_net_prime,a_prime,c_guess,pars)

        call util_c(c_guess,u_c,pars%A_par)
        rho_prime = reshape(u_c(:,1)/u_c(:,2),[M_par,1])


        !We've got plans for most of the variables which we save. The current-period
        !shock index is sim_series%s_ind_sim(t_ind). Take that row of the plans and save it in
        !c_sim_series.
        sim_result%a_prime_sim(t_ind,:) = a_prime(sim_result%s_ind_sim(t_ind),:)
        sim_result%rho_prime_sim(t_ind) = rho_prime(sim_result%s_ind_sim(t_ind),1)
        sim_result%t_prime_sim(t_ind) = t_st_pr
        !(we already have generated and saved a series of shocks)


        sim_result%c_sim(t_ind,:) = c_guess(sim_result%s_ind_sim(t_ind),:)
        sim_result%l_sim(t_ind,:) = lab_mat(sim_result%s_ind_sim(t_ind),:)
        sim_result%u_c_sim(t_ind,:) = u_c(sim_result%s_ind_sim(t_ind),:)
        sim_result%b_net_prime_sim(t_ind,:) = b_net_prime(sim_result%s_ind_sim(t_ind),:)
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




!Get net asset holdings in all periods.
!The following has to be generalized for I>2
    sim_result%b_net_sim(1,:) = pars%b_init(1,2) - pars%b_init(1,1)
    sim_result%b_net_sim(2:pars%T_sim,:) = sim_result%b_net_prime_sim(1:pars%T_sim-1,:)

!Now use the normalization given in parameter file to set an agents asset holdings
!equal to zero.
select case(pars%b_norm)
    case(0)
        !government's assets normalized to zero
        sim_result%B_gov_sim = 0.0_dp
        sim_result%b_sim(:,1) = - pars%mass(1,2)*sim_result%b_net_sim(:,1) !remember that b_net_sim(:,1) is \tilde b_2
        sim_result%b_sim(:,2) =  pars%mass(1,1)*sim_result%b_net_sim(:,1)

        !get the last element of b_prime_sim
        sim_result%b_prime_sim(pars%T_sim,1) = - pars%mass(1,2)*sim_result%b_net_prime_sim(pars%T_sim,1)
        sim_result%b_prime_sim(pars%T_sim,2) = pars%mass(1,1)*sim_result%b_net_prime_sim(pars%T_sim,1)
    case(1)
        !agent 1's assets are 0
        sim_result%b_sim(:,1) = 0.0_dp
        sim_result%b_sim(:,2) =  sim_result%b_net_sim(:,1)
        sim_result%B_gov_sim(:)  = sim_result%b_net_sim(:,1) * pars%mass(1,2)

        !get the last element of b_prime_sim
        sim_result%b_prime_sim(pars%T_sim,1) = 0.0_dp
        sim_result%b_prime_sim(pars%T_sim,2) = sim_result%b_net_sim(pars%T_sim,1)
    case(2)
        !agent 2's assets are 0
        sim_result%b_sim(:,2) = 0.0_dp
        sim_result%b_sim(:,1) =  - sim_result%b_net_sim(:,1)
        sim_result%B_gov_sim(:)  =  - sim_result%b_net_sim(:,1) * pars%mass(1,1)

        !get the last element of b_prime_sim
        sim_result%b_prime_sim(pars%T_sim,1) = - sim_result%b_net_sim(pars%T_sim,1)
        sim_result%b_prime_sim(pars%T_sim,2) = - sim_result%b_net_sim(pars%T_sim,1) * pars%mass(1,1)

    case default
        write(*,*) 'Error: wrong value of b_norm'
        error stop
end select

!get b_prime_sim (just so it's convenient to compute things like transfers and asset income later).
sim_result%b_prime_sim(1:pars%T_sim-1,:) = sim_result%b_sim(2:pars%T_sim,:)


!Get t_sim ('time')
sim_result%t_sim(1) = 0.0_dp !'time' variable in first period is 0
sim_result%t_sim(2:pars%T_sim) = sim_result%t_prime_sim(1:pars%T_sim-1)


    do t_ind=1,pars%T_sim
        !Get t_sim from t_prime_sim

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

        sim_result%trans_to_GDP_sim(t_ind) = sim_result%trans_sim(t_ind) / sim_result%Y_sim(t_ind)

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


        !Compute transfers (from budget constraint of agent 1)
        sim_result%trans_sim(t_ind) = sim_result%c_sim(t_ind,1) + sim_result%b_prime_sim(t_ind,1) - &
        (1.0_dp-sim_result%tau_sim(t_ind))*sim_result%theta_lvl_sim(t_ind,1)*sim_result%l_sim(t_ind,1)
        !the last term depends on R_{t-1}
       if(t_ind > 1) then
            sim_result%trans_sim(t_ind) = sim_result%trans_sim(t_ind) - sim_result%b_sim(t_ind,1)*sim_result%R_sim(t_ind-1)
       else
            sim_result%trans_sim(t_ind) = sim_result%trans_sim(t_ind) - sim_result%b_sim(t_ind,1)
       end if

    end do

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

    !Do the same for transfers (these are the most intersting parts of solution - the others
    !can be computed from the saved average series).
    sim_res_all%trans_first = sim_avg%trans_sim(1)
    sim_res_all%trans_last = sim_avg%trans_sim(pars%T_sim)
    sim_res_all%trans_lastminfirst = sim_res_all%trans_last - sim_res_all%trans_first


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


    !get the quantiles for tau
    do t_ind = 1,pars%T_sim
        !cycle over samples and save all values of a_prime into vector
        do n_ind = 1,pars%N_sim
            tmp_vector(n_ind) = sim_res_all%SRA(n_ind)%tau_sim(t_ind)
        end do

        N = pars%N_sim !number of points in the vector
        nq = 2 !number of quantiles
        q = [0.05_dp, 0.95_dp] !quantile values
        qv = 0.0_dp
        ifail = -1

        call G01AMF(N,tmp_vector,NQ,Q,QV,IFAIL)

        sim_res_all%tau_05(t_ind) = qv(1)
        sim_res_all%tau_95(t_ind) = qv(2)
    end do

    !get the quantiles for transfers
    do t_ind = 1,pars%T_sim
        !cycle over samples and save all values of a_prime into vector
        do n_ind = 1,pars%N_sim
            tmp_vector(n_ind) = sim_res_all%SRA(n_ind)%trans_sim(t_ind)
        end do

        N = pars%N_sim !number of points in the vector
        nq = 2 !number of quantiles
        q = [0.05_dp, 0.95_dp] !quantile values
        qv = 0.0_dp
        ifail = -1

        call G01AMF(N,tmp_vector,NQ,Q,QV,IFAIL)

        sim_res_all%trans_05(t_ind) = qv(1)
        sim_res_all%trans_95(t_ind) = qv(2)
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

    file_path = 'results/'//trim(folder_name)//'/sim/trans.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%trans_sim
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


    file_path = 'results/'//trim(folder_name)//'/sim/b_net.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%b_net_sim
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

    real(dp), dimension(pars%N_sim,pars%T_sim) :: tmp_array
    integer :: sim_ind

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

    file_path = 'results/'//trim(folder_name)//'/sim/stat_tau_05.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%tau_05
    close(unit = 20)


    file_path = 'results/'//trim(folder_name)//'/sim/stat_tau_95.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%tau_95
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_trans_05.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%trans_05
    close(unit = 20)


    file_path = 'results/'//trim(folder_name)//'/sim/stat_trans_95.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%trans_95
    close(unit = 20)

    !save the tax rate for all series
    do sim_ind = 1,pars%N_sim
        tmp_array(sim_ind,:) = sim_res_all%SRA(sim_ind)%tau_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/tau_all.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array
    close(unit = 20)


    !save the gov debt for all series - right now its zeroes but it could become useful later
    !if we for example fix the tax rate...
    do sim_ind = 1,pars%N_sim
        tmp_array(sim_ind,:) = sim_res_all%SRA(sim_ind)%B_gov_to_GDP_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/B_gov_to_GDP_all.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array
    close(unit = 20)

    !save the gov debt for all series - right now its zeroes but it could become useful later
    !if we for example fix the tax rate...
    do sim_ind = 1,pars%N_sim
        tmp_array(sim_ind,:) = sim_res_all%SRA(sim_ind)%B_gov_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/B_gov_all.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array
    close(unit = 20)

    !save transfers for all series
    do sim_ind = 1,pars%N_sim
        tmp_array(sim_ind,:) = sim_res_all%SRA(sim_ind)%trans_to_GDP_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/trans_to_GDP_all.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array
    close(unit = 20)


end subroutine save_stat


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

end module mod_taxhetgr

