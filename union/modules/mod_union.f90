module mod_union
!This module contains type definitions, subroutines, etc. which are written specifically for this project.
!Objects which are useful in other projects are mostly defined in other modules (but these are also different
!between projects).
!(The division is not quite clear-cut, it is just a rough guideline)

use mod_types
use mod_par
use mod_parameters_hard
use mod_utility
use mod_tools

!The following module is for 1-6D bspline.
use bspline_sub_module

implicit none

!type grids_type contains grids. It is written for 2 agents only because
!the addition of third player would make the problem computationally infeasible
!using currently available hardware (going from 5 state variables to 7)
type grids_type
    real(dp),allocatable,dimension(:) :: A1_gr !asset holdings of agent 1 (MU adjusted)
    real(dp),allocatable,dimension(:) :: A2_gr !asset holdings of agent 2 (MU adjusted)
    real(dp),allocatable,dimension(:) :: rho_gr !u1/u2
    real(dp),allocatable,dimension(:) :: aex_gr !MU-adjusted nex foreign debt of country 1.
end type grids_type

!type FP_result_type contains all sorts of results of the first-period problem
type FP_result_type
    real(dp), dimension(1,I_par) :: c,l,b_prime,a_prime,u_c,u_l
    real(dp) :: rho_prime,aex_prime,b_ex_prime
    integer :: s_ind_init !initial index - this is usually drawn from a stationary distribution
    !(so that we don't get a solution dependent on initial state which can imply large quantitative
    !trend in government expenditure in expectation).
end type FP_result_type

!Type sim_result type contains
type sim_result_type
    !First index is time period, 2nd is agent
    real(dp), allocatable, dimension(:,:) :: c_sim,l_sim,b_sim,b_prime_sim,a_prime_sim,u_c_sim,u_l_sim,tau_sim,U_sim
    !Note: there is no a_prime and rho_prime in the first period.
    real(dp), allocatable, dimension(:) :: rho_prime_sim,aex_prime_sim,b_ex_sim,b_ex_prime_sim,B_ex_to_GDP_sim
    real(dp), allocatable, dimension(:) :: R_sim,SWF_sim,term_sim,totret_sim
    real(dp), allocatable, dimension(:) :: loss_sim !can be useful for debugging.

    integer, allocatable, dimension(:) :: s_ind_sim

    real(dp), allocatable, dimension(:,:) :: G_sim,G_to_Y_sim
    real(dp), allocatable, dimension(:,:) :: Y_sim,B_to_GDP_sim


    real(dp) :: a_slack_all_avg !average 'slackness' - basically low values (particularly negative) mean
    !that we have an issue in the sense that the a_prime state is close to bounds often, or even
    !gets outside of the bounds. The lower this is the worse the solution is expected to be and
    !this can then be used to eliminate 'outliers'.

end type sim_result_type

!Type sim_result_all type contains all N_sim simulated series, and some addtional statistics computed
!locally (other statistics will be computed by Matlab because it is more convenient to code it).
type sim_result_all_type
    type(sim_result_type), dimension(:), allocatable :: SRA !sim_result_all - contains all simulated series.

    !average production in each country for each shock realization - this is useful for calibration, such as when
    !we want to set the magnitude of shocks such that production drops by a certain percentage in
    !a recession (when g low)
    real(dp), dimension(M_par,I_par) :: avg_Y

    !Also save the 0.05 and 0.95 quantile of realizations of state variables. This is important for detecting
    !issues where the bounds of grid are wrongly chosen (the mean could still be within the bounds but
    !for some particular simulations (shock realizations the bounds could be hit or we could
    !be getting really close to them.
    !These definitions assume I=2, change later if needed!!!
    real(dp), dimension(:), allocatable :: rho_prime_05, rho_prime_95, aex_prime_05,aex_prime_95
    real(dp), dimension(:,:), allocatable :: a_prime_05, a_prime_95
end type sim_result_all_type

!Type SS_type contains a deterministic steady-state allocation
type SS_type
    real(dp), dimension(1,I_par) :: c= 0.0_dp ,l = 0.0_dp !steady-state consumption and labour supply (one per country)
    real(dp), dimension(1,I_par) :: util = 0.0_dp !steady-state per-period utility (in each country)
    real(dp), dimension(1,I_par) :: b= 0.0_dp ,a =0.0_dp !steady-state asset holdings (in levels and MU-adjusted)
    real(dp) :: b_ex= 0.0_dp ,a_ex= 0.0_dp  !steady-state net external debt of country 1
    real(dp) :: rho = 0.0_dp !steady-state ratio of marginal utilities

    real(dp) :: c_ss_tot !total consumption
    real(dp) :: kappa_g2ming1 !the excess of government expenditures in country 2 over country 1 in the SS,
    !multiplied by LFFC_kappa. This is used when improving a naive initial guess in LFFC).
end type SS_type


!Shared data in sim_series (we will declare a coarray of this type in the main program)
type sim_shared_type
    real(dp), dimension(I_par) :: a_st
    real(dp), dimension(I_par-1) :: rho_st
    real(dp)  :: aex_st
    integer :: g_st_ind
    !choice variables - initial guesses
    real(dp), dimension(ncvar_par) :: xc_last,xc_cent,xc_interp,xc_cent_loc,xc
    !bset choices found by an image - really I don't need an extra copy but it's easier this way.
    real(dp) :: bcsf_fc,fc
    real(dp), dimension(ncvar_par) :: bcsf_x

    real(dp), dimension(5,ncvar_par) :: xc_cent_all !this is an array of choices around which we generate the permutations.
    !so far I assume that there are at most 5 different choices around which we permutate. Increase this later if needed.
    integer :: xc_cent_count !the first one is the number of choices that we use. The second one is an index
    !for cycling over these choices.
end type sim_shared_type

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
    sim_plus_sim%a_prime_sim = sr1%a_prime_sim + sr2%a_prime_sim
    sim_plus_sim%u_c_sim = sr1%u_c_sim + sr2%u_c_sim
    sim_plus_sim%u_l_sim = sr1%u_l_sim + sr2%u_l_sim
    sim_plus_sim%tau_sim = sr1%tau_sim + sr2%tau_sim
    sim_plus_sim%rho_prime_sim = sr1%rho_prime_sim + sr2%rho_prime_sim
    sim_plus_sim%aex_prime_sim = sr1%aex_prime_sim + sr2%aex_prime_sim
    sim_plus_sim%b_ex_sim = sr1%b_ex_sim + sr2%b_ex_sim
    sim_plus_sim%b_ex_prime_sim = sr1%b_ex_prime_sim + sr2%b_ex_prime_sim
    sim_plus_sim%B_ex_to_GDP_sim = sr1%B_ex_to_GDP_sim + sr2%B_ex_to_GDP_sim
    sim_plus_sim%loss_sim = sr1%loss_sim + sr2%loss_sim
    sim_plus_sim%s_ind_sim = sr1%s_ind_sim + sr2%s_ind_sim
    sim_plus_sim%G_sim = sr1%G_sim + sr2%G_sim
    sim_plus_sim%G_to_Y_sim = sr1%G_to_Y_sim + sr2%G_to_Y_sim
    sim_plus_sim%Y_sim = sr1%Y_sim + sr2%Y_sim
    sim_plus_sim%B_to_GDP_sim = sr1%B_to_GDP_sim + sr2%B_to_GDP_sim
    sim_plus_sim%a_slack_all_avg = sr1%a_slack_all_avg + sr2%a_slack_all_avg
    sim_plus_sim%R_sim = sr1%R_sim + sr2%R_sim
    sim_plus_sim%u_sim = sr1%u_sim + sr2%u_sim
    sim_plus_sim%swf_sim = sr1%swf_sim + sr2%swf_sim
    sim_plus_sim%term_sim = sr1%term_sim + sr2%term_sim
    sim_plus_sim%totret_sim = sr1%totret_sim + sr2%totret_sim


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
    sim_min_sim%a_prime_sim = sr1%a_prime_sim - sr2%a_prime_sim
    sim_min_sim%u_c_sim = sr1%u_c_sim - sr2%u_c_sim
    sim_min_sim%u_l_sim = sr1%u_l_sim - sr2%u_l_sim
    sim_min_sim%tau_sim = sr1%tau_sim - sr2%tau_sim
    sim_min_sim%rho_prime_sim = sr1%rho_prime_sim - sr2%rho_prime_sim
    sim_min_sim%aex_prime_sim = sr1%aex_prime_sim - sr2%aex_prime_sim
    sim_min_sim%b_ex_sim = sr1%b_ex_sim - sr2%b_ex_sim
    sim_min_sim%b_ex_prime_sim = sr1%b_ex_prime_sim - sr2%b_ex_prime_sim
    sim_min_sim%B_ex_to_GDP_sim = sr1%B_ex_to_GDP_sim - sr2%B_ex_to_GDP_sim
    sim_min_sim%loss_sim = sr1%loss_sim - sr2%loss_sim
    sim_min_sim%s_ind_sim = sr1%s_ind_sim - sr2%s_ind_sim !Ok, this is defined, but it's useless to compute averages for this
    sim_min_sim%G_sim = sr1%G_sim - sr2%G_sim
    sim_min_sim%G_to_Y_sim = sr1%G_to_Y_sim - sr2%G_to_Y_sim
    sim_min_sim%Y_sim = sr1%Y_sim - sr2%Y_sim
    sim_min_sim%B_to_GDP_sim = sr1%B_to_GDP_sim - sr2%B_to_GDP_sim
    sim_min_sim%a_slack_all_avg = sr1%a_slack_all_avg - sr2%a_slack_all_avg
    sim_min_sim%R_sim = sr1%R_sim - sr2%R_sim
    sim_min_sim%u_sim = sr1%u_sim - sr2%u_sim
    sim_min_sim%swf_sim = sr1%swf_sim - sr2%swf_sim
    sim_min_sim%term_sim = sr1%term_sim - sr2%term_sim
    sim_min_sim%totret_sim = sr1%totret_sim - sr2%totret_sim

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
    sim_multi%a_prime_sim = sr%a_prime_sim * multi
    sim_multi%u_c_sim = sr%u_c_sim * multi
    sim_multi%u_l_sim = sr%u_l_sim * multi
    sim_multi%tau_sim = sr%tau_sim * multi
    sim_multi%rho_prime_sim = sr%rho_prime_sim * multi
    sim_multi%aex_prime_sim = sr%aex_prime_sim * multi
    sim_multi%b_ex_sim = sr%b_ex_sim * multi
    sim_multi%b_ex_prime_sim = sr%b_ex_prime_sim * multi
    sim_multi%B_ex_to_GDP_sim = sr%B_ex_to_GDP_sim * multi
    sim_multi%loss_sim = sr%loss_sim * multi
    sim_multi%s_ind_sim = sr%s_ind_sim * multi !Ok, this is defined, but it's useless to compute averages for this
    sim_multi%G_sim = sr%G_sim * multi
    sim_multi%G_to_Y_sim = sr%G_to_Y_sim * multi
    sim_multi%Y_sim = sr%Y_sim * multi
    sim_multi%B_to_GDP_sim = sr%B_to_GDP_sim * multi
    sim_multi%a_slack_all_avg = sr%a_slack_all_avg * multi
    sim_multi%R_sim = sr%R_sim * multi
    sim_multi%u_sim = sr%u_sim * multi
    sim_multi%swf_sim = sr%swf_sim * multi
    sim_multi%term_sim = sr%term_sim * multi
    sim_multi%totret_sim = sr%totret_sim * multi
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

    div_inverse = 1.0_dp/div

    sim_div = sim_multi(sr,div_inverse)

end function sim_div

!Subroutine b_init_perm generates a random permutation of the initial conditions.
subroutine b_init_permutate(b_init,b_init_perm,b_ex_init,b_ex_init_perm,perm_eps)
    real(dp), dimension(1,I_par), intent(in) :: b_init !initial conditions as given in parameter file
    real(dp), intent(in) :: b_ex_init
    real(dp), dimension(1,I_par), intent(out) :: b_init_perm !permutated initial conditions
    real(dp), intent(out) :: b_ex_init_perm

    real(dp) :: perm_eps

    real(dp), dimension(3) :: rnd !random number drawn from 0,1 interval (uniform distribution)

    call random_number(rnd)
    b_init_perm(1,1) = b_init(1,1) + perm_eps*(2.0_dp*rnd(1) - 1.0_dp)

    b_init_perm(1,2) = b_init(1,2) + perm_eps*(2.0_dp*rnd(2) - 1.0_dp)

    b_ex_init_perm = b_ex_init + perm_eps*(2.0_dp*rnd(3) - 1.0_dp)

end subroutine b_init_permutate


!Subroutine c_inc_mtv transforms the incomplete consumption matrix (with element
!M_par,I_par missing) into a row form, so it can be used by some NAG functions
!for maximization. It just reshapes the matrix (using standard Fortran ordering
!where first index varies most rapidly.
subroutine c_inc_mtv(c_inc_mat,c_inc_vec)
    real(dp),dimension(M_par,I_par), intent(in) :: c_inc_mat
    real(dp), dimension(M_par*I_par - 1), intent(out) :: c_inc_vec
    real(dp), dimension(M_par*I_par) :: c_inc_vec_tmp !one more element

    c_inc_vec_tmp = reshape(c_inc_mat,[M_par*I_par])
    c_inc_vec = c_inc_vec_tmp(1:M_par*I_par - 1)

    !ncvar_par = M_par*I_par - 1 is the number of choice variables
end subroutine c_inc_mtv


!Subroutine check_constr checks whether constraints of the problem are violated, and it returns a value
!of loss function (which can potentially be adjusted downwards if some constraints are strictly satisfied -
!this is used in finding feasible choice to increase the chance that we are not at the edge of a feasible
!set when starting maximization).
subroutine check_constr(pars,c,l,a_pr,rho_pr,aex_pr,constr_ok,loss,rofv,neg_loss)
    type(par), intent(in) :: pars
    real(dp), dimension(M_par,I_par), intent(in) :: c,l,a_pr !consumption, labour, next-period debt
    real(dp), dimension(M_par,1), intent(in) :: rho_pr, aex_pr !next period states rho and aex

    logical, intent(in) :: neg_loss !If this is true then loss is negative for points which strictly satisfy
    !all constraints, and it is 'more negative' the further from bounds the point is. This is used to find
    !the 'most feasible' point in the sense of laying in the middle of some large convex feasible set (we need
    !to do this because of nonconvexity of feasible sets).

    logical, intent(out) :: constr_ok !this will be true if all constraints are satisfied
    real(dp),intent(out) :: loss !loss stemming from violating constraints.

    logical, intent(in) :: rofv !return on first violation
    !If this is true, once a constraint is violated, the control is returned, and
    !we don't finish checking the other constraints. This should be false when looking for
    !feasible choices until we find one feasible choice, then it should be true
    !(no point looking for not strictly feasible choices if we already have a strictly
    !feasible one - that is the first preference). But this is the case only if we are trying
    !random points. If we are using an improvement algorithm (minimizing loss), we still need to
    !compute violation of all constraints to get the loss. Need to be careful when using value true,
    !but in some cases we can save a lot of time due to this.

    !Locally used variables.
    real(dp), parameter, dimension(M_par,1) :: ones = 1.0_dp

    !Parameters for loss function computation are given in pars.

    !initialize constraints_ok and loss
    constr_ok = .true.
    loss = 0.0_dp

    !(1) Check whether consumption or labour supply are negative. If this is the case than there is an error in the program,
    !because this should be checked for earlier by subroutine get_av.
    !In that case we do not recover the remaining variables,
    !and we should not call check_constr (rather apply a fixed penalty) or evaluate the current return.
    if(any(c<=0.0_dp) .or. any(l<0.0_dp)) then
        write(*,*) 'Error: negative consumption or labour supply in subroutine check_const. This should never happen!'
        write(*,*) 'c = ',c
        write(*,*) 'l = ',l
        error stop
    end if

    !(2) Check that labour supply does not exceed time endowment. If yes, apply a convex cost.
    if(any(l>pars%l_max_mat)) then
        constr_ok = .false.
        if(rofv) return !return before even computing loss because if rofv is true,
            !we only care about strict feasibility and not about value of loss function.

        !sum of squared violations multiplied by a parameter.
        loss = loss + pars%loss3a*sum(max(l-pars%l_max_mat,0.0_dp)**2)
    end if

    !(3) Checking constraints on a_pr.

    !Here loss function will be particularly important (no jumps)
    if(any(a_pr>pars%a_max_mat)) then
        constr_ok = .false.
        if(rofv) return
        loss = loss + pars%loss1a*sum(max(a_pr-pars%a_max_mat,0.0_dp)**2)
    end if
    if(any(a_pr<pars%a_min_mat)) then
        constr_ok = .false.
        if(rofv) return
        loss = loss + pars%loss1a*sum(max(pars%a_min_mat-a_pr,0.0_dp)**2)
    end if

    !(4) Checking constraints on rho_pr
    if(any(rho_pr>pars%rho_max_vec)) then
        constr_ok = .false.
        if(rofv) return
        loss = loss + pars%loss2a*sum(max(rho_pr-pars%rho_max_vec,0.0_dp)**2)
    end if
    if(any(rho_pr<pars%rho_min_vec)) then
        constr_ok = .false.
        if(rofv) return
        loss = loss + pars%loss2a*sum(max(pars%rho_min_vec-rho_pr,0.0_dp)**2)
    end if

    !(5) Checking constraints on aex_pr
    if(any(aex_pr>pars%aex_max_vec)) then
        constr_ok = .false.
        if(rofv) return
        loss = loss + pars%loss1a*sum(max(aex_pr-pars%aex_max_vec,0.0_dp)**2)
    end if
    if(any(aex_pr<pars%aex_min_vec)) then
        constr_ok = .false.
        if(rofv) return
        loss = loss + pars%loss1a*sum(max(pars%aex_min_vec-aex_pr,0.0_dp)**2)
    end if

    !If the loss is not zero, then increase it by the loss constant. The loss constant
    !is usually zero, but sometimes it is positive to make sure that no 'cheating' by
    !exploiting the fact that states are truncated is advantageous. However, this introduces
    !a discontinuity, so it is probably safer to just use higher values of the coefficients
    !loss1a,loss2a,loss3a. However, leave this  option in the code for testing of implications.
    if(loss > 0.000000000001_dp) then
        loss = loss + pars%loss_const
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

    !For a_prime scale the adjustment by the range. Otherwise the terms here will grow linearly with the width
    !of feasible set for MU adjusted asset holdings. It could then possibly outweigh loss and make us choose
    !points which are in the middle of feasible set for a_prime but are infeasible with respect to
    !other constraints. Also it would favour feasibility in sense of a_prime relative to feasibility in sense of c or l.
    !The normalization assumes that all elements of a_max_mat - a_min_mat are positive. If not we would
    !get all sorts of problems, crashes, and non-sensical results elsewhere.
    loss = loss - sqrt(sum(max((a_pr-pars%a_min_mat)/(pars%a_max_mat-pars%a_min_mat),0.0_dp)**2))
    loss = loss - sqrt(sum(min((a_pr-pars%a_max_mat)/(pars%a_max_mat-pars%a_min_mat),0.0_dp)**2))

    !Do the same for rho_prime. This is relevant only in case of stochastic calibration or
    !corner solutions. In the case of deterministic calibration and interior solutions the value
    !of rho is fixed so it will be just a constant term in objective function of a LFFC loss minimization
    !subroutine.
    loss = loss - sqrt(sum(max((rho_pr-pars%rho_min_vec)/(pars%rho_max_vec-pars%rho_min_vec),0.0_dp)**2))
    loss = loss - sqrt(sum(min((rho_pr-pars%rho_max_vec)/(pars%rho_max_vec-pars%rho_min_vec),0.0_dp)**2))

    !and for aex_prime
    loss = loss - sqrt(sum(max((aex_pr-pars%aex_min_vec)/(pars%aex_max_vec-pars%aex_min_vec),0.0_dp)**2))
    loss = loss - sqrt(sum(min((aex_pr-pars%aex_max_vec)/(pars%aex_max_vec-pars%aex_min_vec),0.0_dp)**2))
end if

    !Cap loss at the upper bound defined in mod_parameters_hard. This is just to avoid bugs. And mitigate divergence
    !issues (there is no point trying local improvements if loss is so huge - it's better to just start with a new guess
    !somewhere else).
    loss = min(loss,max_loss)

end subroutine check_constr


!Subroutine CPC_ret computes (-1)*( current + continuation return), given current state
!and the incomplete consumption matrix (in vector form).
!
!It is written in such form that it can be used by NAG Fortran Subroutine E04JYF
!and possibly other subroutines. For this reason, many of the input parameters
!are passed to this function using common block.
subroutine CPC_ret(n,xc,fc,iuser,ruser)
        real (dp), intent (out) :: fc
        integer, intent (in) :: n
        real (dp), intent (inout) :: ruser(*)
        real (dp), intent (in) :: xc(n)
        integer, intent (inout) :: iuser(*)

        real(dp), dimension(:,:,:,:,:), pointer :: V_old_pntr !pointer to V_old, passed using common block
        real(dp), dimension(:), pointer :: a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr !pointers so all grids

        !pointer pars_pntr points to parameters type (it is more convenient to pass it using a common
        !block than create copies of everything we need and put them in a long common block). We
        !just need to be careful to not make any assignment statements because we would then modify
        !parameters elsewhere in the program!!!
        type(par), pointer :: pars_pntr

        !Values of state variables as saved in the grid (hence the suffix st).
        real(dp), dimension(I_par) :: a_st
        real(dp), dimension(I_par-1) :: rho_st
        real(dp) :: aex_st
        integer :: g_st_ind

        !Allocation and next-period states
        real(dp), dimension(M_par,I_par) :: c,l
        real(dp), dimension(M_par,I_par) :: a_pr
        real(dp), dimension(M_par,1) :: rho_pr, aex_pr

        logical :: constr_ok

        integer :: g_ind !index of shock realization

        real(dp) :: loss !loss from violating constraints.

        real(dp), dimension(M_par,1) :: conval_all !continuation value in all states

        integer :: getavfail

        integer, dimension(M_par,4) :: j_guess !The initial guess of position of next-period states
        !in the grid, for evry shock realization.

        !It is very important that nothing that we pass to this function using common block
        !is changed here! Otherwise it might lead to errors elsewhere in the program.
        common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
        !to some subroutines.
        common /pars_pointer/ pars_pntr
        common /grids/ a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr !grids

        common /cpc_ret_bloc/ a_st,rho_st,aex_st,g_st_ind,j_guess

        common /commloss/ loss !For debugging purposes only (so that loss is reported one level higher, not just fc)

        !For the choice variables xc, get all the remaining variables using subroutine getav
        !Compute all the remaining variables implied by the choice xc
        call get_av(pars_pntr,xc,a_st,rho_st,aex_st,g_st_ind,getavfail,c,l,a_pr,rho_pr,aex_pr,.false.)

        !If non-negativity of consumption or labour is violated, then current return is not defined. Set fc = max_loss
        !(this should guarantee such a point is never deemed optimal) and return.
        if(getavfail > 0) then
            fc = max_loss !maximum loss. This introduces a significant discontinuity, need to be careful about implications elsewhere
            return
        end if

        !Compute the loss (unadjusted!!!)
        call check_constr(pars_pntr,c,l,a_pr,rho_pr,aex_pr,constr_ok,loss,.false.,.false.)

        !Now that we know consumption and labour plan, we can compute the current (expected) return
        fc = SWF(c,l,g_st_ind,pars_pntr%P,pars_pntr%mass,pars_pntr%alpha,pars_pntr%A_par,&
        pars_pntr%B_par,pars_pntr%sigma_par,pars_pntr%Gamma_par)

        !Now we have current return in fc, and loss in loss. We just need to compute the continuation return.

        !If the choice is not strictly feasible, force it to be within bounds (extrapolation
        !is very bad, usually leads to divergence). So effectively if the algorithm chooses a point
        !outside of the grid boundaries, we use neares neighbour extapolation (sort of) and add
        !a loss. This is much more stable (otherwise in some corners of the grid the loss if we
        !use extrapolation is compounded and the program diverges). This puts a lower bound on it and it is
        !fine as long as the region where this happens was not optimally reached in the first place (need
        !to experiment with loss).
        if(.not. constr_ok) then
            a_pr=max(pars_pntr%a_min_mat,a_pr)
            a_pr=min(pars_pntr%a_max_mat,a_pr)
            rho_pr=max(pars_pntr%rho_min_vec,rho_pr)
            rho_pr=min(pars_pntr%rho_max_vec,rho_pr)
            aex_pr=max(pars_pntr%aex_min_vec,aex_pr)
            aex_pr=min(pars_pntr%aex_max_vec,aex_pr)
        end if


        !Once we have the point forced within the grid boundaries, get the continuation
        !value by interpolation. We need to get a separate continuation value for all
        !states.

        !Which interpolation function we use depends on the value of VFI_interp_mode (set in parameter file).
        select case(pars_pntr%VFI_interpolation_mode)
            case(1) !Quadrilinear interpolation
                do g_ind = 1,M_par
                    call interp_V(V_old_pntr,[a_pr(g_ind,1),a_pr(g_ind,2),rho_pr(g_ind,1),aex_pr(g_ind,1)],j_guess(g_ind,:)&
                    ,min(g_ind,pars_pntr%M_grid),conval_all(g_ind,1),a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr,size(a1_gr_pntr),&
                    size(rho_gr_pntr),size(aex_gr_pntr))
                end do
            case(2) !Local Shepard's interpolation
                do g_ind = 1,M_par
                    call interp_V_shep(V_old_pntr,[a_pr(g_ind,1),a_pr(g_ind,2),rho_pr(g_ind,1),aex_pr(g_ind,1)],j_guess(g_ind,:)&
                    ,min(g_ind,pars_pntr%M_grid),conval_all(g_ind,1),a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr,size(a1_gr_pntr),&
                    size(rho_gr_pntr),size(aex_gr_pntr),pars_pntr%shep_n,pars_pntr%shep_norm,pars_pntr%shep_norm_val,pars_pntr%shep_floor)
                end do
            case(3) !bspline
                do g_ind = 1,M_par
                    call interp_V_spline(V_old_pntr,[a_pr(g_ind,1),a_pr(g_ind,2),rho_pr(g_ind,1),aex_pr(g_ind,1)],&
                    min(g_ind,pars_pntr%M_grid),conval_all(g_ind,1),a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr,size(a1_gr_pntr),&
                    size(rho_gr_pntr),size(aex_gr_pntr),pars_pntr%bspline_k_a1,pars_pntr%bspline_k_a2,pars_pntr%bspline_k_rho,pars_pntr%bspline_k_aex)
                end do
        end select

        !Now we have continuation value for every possible shock realization. Add the discounted
        !expectation of this term to the current return. (g_st_ind is the last shock realization)

        fc = fc + pars_pntr%beta*maxval(matmul(pars_pntr%P(g_st_ind,:),conval_all))
        !maxval is there simply to convert (1,1) array to scalar.

        !Now fc continues (expected) current plus continuation return for the choice xc

        !Subtract loss stemming from violating some constraints, and multiply by -1 (because we are using a minimization
        !subroutine)
        fc = (fc - loss)*(-1.0_dp)
end subroutine CPC_ret

!This is a copy of CPC_ret but with advanced interface (iuser is used),
!and more functions. This subroutine should currenly be used only for debugging purposes,
!and to compute terminal values in simulation (not in VFI).
!
!The subroutine does the following depending on iuser:
!iuser(1) determines what shall be returned in the scalar fc
!
!it never makes sense to look for a choice which minimizes the current return,
!   0: CPC_ret behaviour (current plus continuation return including all penalizations times -1)
!   1: current return only (times -1)
!   2: continuation return only (times -1). Discounted.
!   3: loss
!   4: distance penalization (this is zero unless iuser(3) > 0.
!   (note: the penalizations are not returned as -1 multiple because we sometimes might want to see what
!   choice minimizes them. For all the other returnable options we want to maximize them always if any
!   optimization is made there).
!iuser(2) determines which state the variables returned correpond to
!   0: expectation
!   s > 0: shock realization s
!   (note: this only affects the current return (current-period SWF) and continuation value. Loss
!   and distance penalization are not computed for individual states).
!iuser(3) determines whether penalization should be done for distance of next-period states from a given
!   point in state space.
!   0: no penalization
!   1: penalization for distance from last-period state
!   2: penaliation for distance from the center of the grid
!   (note: These variables can be used for debugging purposes only such as figuring out whether it is
!   even possible to get to the 'middle of the grid' from some parts of the state space. If used in actual
!   simulations or VFI to compute the optimal policy, this will invalidate the results).
!In all cases - if get_av_fail, then the function returns the value max_loss in fc.
subroutine CPC_ret_adv(n,xc,fc,iuser,ruser)
        real (dp), intent (out) :: fc
        integer, intent (in) :: n
        real (dp), intent (inout) :: ruser(*)
        real (dp), intent (in) :: xc(n)
        integer, intent (inout) :: iuser(*)

        real(dp), dimension(:,:,:,:,:), pointer :: V_old_pntr !pointer to V_old, passed using common block
        real(dp), dimension(:), pointer :: a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr !pointers so all grids

        !pointer pars_pntr points to parameters type (it is more convenient to pass it using a common
        !block than create copies of everything we need and put them in a long common block). We
        !just need to be careful to not make any assignment statements because we would then modify
        !parameters elsewhere in the program!!!
        type(par), pointer :: pars_pntr

        !Values of state variables as saved in the grid (hence the suffix st).
        real(dp), dimension(I_par) :: a_st
        real(dp), dimension(I_par-1) :: rho_st
        real(dp) :: aex_st
        integer :: g_st_ind

        !Allocation and next-period states
        real(dp), dimension(M_par,I_par) :: c,l
        real(dp), dimension(M_par,I_par) :: a_pr
        real(dp), dimension(M_par,1) :: rho_pr, aex_pr

        logical :: constr_ok

        integer :: g_ind !index of shock realization

        real(dp) :: loss !loss from violating constraints.

        real(dp), dimension(M_par,1) :: conval_all !continuation value in all states

        integer :: getavfail

        integer, dimension(M_par,4) :: j_guess !The initial guess of position of next-period states
        !in the grid, for evry shock realization.

        !It is very important that nothing that we pass to this function using common block
        !is changed here! Otherwise it might lead to errors elsewhere in the program.
        common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
        !to some subroutines.
        common /pars_pointer/ pars_pntr
        common /grids/ a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr !grids

        common /cpc_ret_bloc/ a_st,rho_st,aex_st,g_st_ind,j_guess

        common /commloss/ loss !For debugging purposes only (so that loss is reported one level higher, not just fc)

        if(this_image() == 1) then
            if(iuser(1) /= 2) then
                write(*,*) 'Error in subroutine CPC_ret_adv. At this point it only works'
                write(*,*) 'if iuser(1) == 2 - which computes the continuation value.'
            end if
        end if

        !This subroutine is not written with efficiency in mind. Basically it runs just like CPC_ret
        !but at some points I compute extra information and I possibly return (overwriting the current
        !value of fc).

        !For the choice variables xc, get all the remaining variables using subroutine getav
        !Compute all the remaining variables implied by the choice xc
        call get_av(pars_pntr,xc,a_st,rho_st,aex_st,g_st_ind,getavfail,c,l,a_pr,rho_pr,aex_pr,.false.)

        !If I want to compute current return only(iuser(1) == 1) I need to call util
        !and manually compute the value for each state, and then return one of them.
        !esentially just like what I do below for the continuation return.

        !If non-negativity of consumption or labour is violated, then current return is not defined. Set fc = max_loss
        !(this should guarantee such a point is never deemed optimal) and return.
        if(getavfail > 0) then
            fc = max_loss !maximum loss. This introduces a significant discontinuity, need to be careful about implications elsewhere
            return
        end if

        !Compute the loss (unadjusted!!!)
        call check_constr(pars_pntr,c,l,a_pr,rho_pr,aex_pr,constr_ok,loss,.false.,.false.)

        !Now that we know consumption and labour plan, we can compute the current (expected) return
        fc = SWF(c,l,g_st_ind,pars_pntr%P,pars_pntr%mass,pars_pntr%alpha,pars_pntr%A_par,&
        pars_pntr%B_par,pars_pntr%sigma_par,pars_pntr%Gamma_par)

        !Now we have current return in fc, and loss in loss. We just need to compute the continuation return.

        !If the choice is not strictly feasible, force it to be within bounds (extrapolation
        !is very bad, usually leads to divergence). So effectively if the algorithm chooses a point
        !outside of the grid boundaries, we use neares neighbour extapolation (sort of) and add
        !a loss. This is much more stable (otherwise in some corners of the grid the loss if we
        !use extrapolation is compounded and the program diverges). This puts a lower bound on it and it is
        !fine as long as the region where this happens was not optimally reached in the first place (need
        !to experiment with loss).
        if(.not. constr_ok) then
            a_pr=max(pars_pntr%a_min_mat,a_pr)
            a_pr=min(pars_pntr%a_max_mat,a_pr)
            rho_pr=max(pars_pntr%rho_min_vec,rho_pr)
            rho_pr=min(pars_pntr%rho_max_vec,rho_pr)
            aex_pr=max(pars_pntr%aex_min_vec,aex_pr)
            aex_pr=min(pars_pntr%aex_max_vec,aex_pr)
        end if


        !Once we have the point forced within the grid boundaries, get the continuation
        !value by interpolation. We need to get a separate continuation value for all
        !states.

        !Which interpolation function we use depends on the value of VFI_interp_mode (set in parameter file).
        select case(pars_pntr%VFI_interpolation_mode)
            case(1) !Quadrilinear interpolation
                do g_ind = 1,M_par
                    call interp_V(V_old_pntr,[a_pr(g_ind,1),a_pr(g_ind,2),rho_pr(g_ind,1),aex_pr(g_ind,1)],j_guess(g_ind,:)&
                    ,min(g_ind,pars_pntr%M_grid),conval_all(g_ind,1),a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr,size(a1_gr_pntr),&
                    size(rho_gr_pntr),size(aex_gr_pntr))
                end do
            case(2) !Local Shepard's interpolation
                do g_ind = 1,M_par
                    call interp_V_shep(V_old_pntr,[a_pr(g_ind,1),a_pr(g_ind,2),rho_pr(g_ind,1),aex_pr(g_ind,1)],j_guess(g_ind,:)&
                    ,min(g_ind,pars_pntr%M_grid),conval_all(g_ind,1),a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr,size(a1_gr_pntr),&
                    size(rho_gr_pntr),size(aex_gr_pntr),pars_pntr%shep_n,pars_pntr%shep_norm,pars_pntr%shep_norm_val,pars_pntr%shep_floor)
                end do
            case(3) !bspline
                do g_ind = 1,M_par
                    call interp_V_spline(V_old_pntr,[a_pr(g_ind,1),a_pr(g_ind,2),rho_pr(g_ind,1),aex_pr(g_ind,1)],&
                    min(g_ind,pars_pntr%M_grid),conval_all(g_ind,1),a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr,size(a1_gr_pntr),&
                    size(rho_gr_pntr),size(aex_gr_pntr),pars_pntr%bspline_k_a1,pars_pntr%bspline_k_a2,pars_pntr%bspline_k_rho,pars_pntr%bspline_k_aex)
                end do
        end select

        if(iuser(1) == 2) then
            !In this case we want to compute and return the continuation value
            if(iuser(2) == 0) then
                !we compute the expected continuation return (discounted)
                fc = -pars_pntr%beta*maxval(matmul(pars_pntr%P(g_st_ind,:),conval_all))
                return
            else
                fc = -pars_pntr%beta * conval_all(iuser(2),1)
                return
            end if
        end if

        !Now we have continuation value for every possible shock realization. Add the discounted
        !expectation of this term to the current return. (g_st_ind is the last shock realization)

        fc = fc + pars_pntr%beta*maxval(matmul(pars_pntr%P(g_st_ind,:),conval_all))
        !maxval is there simply to convert (1,1) array to scalar.

        !Now fc continues (expected) current plus continuation return for the choice xc

        !Subtract loss stemming from violating some constraints, and multiply by -1 (because we are using a minimization
        !subroutine)
        fc = (fc - loss)*(-1.0_dp)
end subroutine CPC_ret_adv


!Subroutine CPC_ret2 is a copy of CPC_ret, the only difference being the presence of one
!more output (inform). This is defined for use in a different NAg subroutine than CPC_ret.
subroutine CPC_ret2(n,xc,fc,iuser,ruser,inform)
        real (dp), intent (out) :: fc
        integer, intent (in) :: n
        real (dp), intent (inout) :: ruser(*)
        real (dp), intent (in) :: xc(n)
        integer, intent (inout) :: iuser(*) !In this instance, it is used to pass what_binds
        integer, intent(out) :: inform

        real(dp), dimension(:,:,:,:,:), pointer :: V_old_pntr !pointer to V_old, passed using common block
        real(dp), dimension(:), pointer :: a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr !pointers so all grids

        !pointer pars_pntr points to parameters type (it is more convenient to pass it using a common
        !block than create copies of everything we need and put them in a long common block). We
        !just need to be careful to not make any assignment statements because we would then modify
        !parameters elsewhere in the program!!!
        type(par), pointer :: pars_pntr

        !Values of state variables as saved in the grid (hence the suffix st).
        real(dp), dimension(I_par) :: a_st
        real(dp), dimension(I_par-1) :: rho_st
        real(dp) :: aex_st
        integer :: g_st_ind

        !Allocation and next-period states
        real(dp), dimension(M_par,I_par) :: c,l
        real(dp), dimension(M_par,I_par) :: a_pr
        real(dp), dimension(M_par,1) :: rho_pr, aex_pr

        logical :: constr_ok

        integer :: g_ind !index of shock realization

        real(dp) :: loss !loss from violating constraints.

        real(dp), dimension(M_par,1) :: conval_all !continuation value in all states

        integer :: getavfail

        integer, dimension(M_par,4) :: j_guess !The initial guess of position of next-period states
        !in the grid, for evry shock realization.

        !It is very important that nothing that we pass to this function using common block
        !is changed here! Otherwise it might lead to errors elsewhere in the program.
        common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
        !to some subroutines.
        common /pars_pointer/ pars_pntr
        common /grids/ a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr !grids

        common /cpc_ret_bloc/ a_st,rho_st,aex_st,g_st_ind,j_guess

        common /commloss/ loss !For debugging purposes only (so that loss is reported one level higher, not just fc)

        inform = 0 !initialize inform. If this is returned as a negative number, then the optimization
        !subroutine stops. Do not use this, leads to bugs. If this is used we would need careful handling of errors
        !one level above.

        !For the choice variables xc, get all the remaining variables using subroutine getav
        !Compute all the remaining variables implied by the choice xc
        call get_av(pars_pntr,xc,a_st,rho_st,aex_st,g_st_ind,getavfail,c,l,a_pr,rho_pr,aex_pr,.false.)

        !If non-negativity of consumption or labour is violated, then current return is not defined. Set fc = max_loss
        !(this should guarantee such a point is never deemed optimal) and return.
        if(getavfail > 0) then
            fc = max_loss !maximum loss. This introduces a significant discontinuity, need to be careful about implications elsewhere
            !inform = -1 DO NOT USE THIS! If this is used, then on exit, sometimed fmax does not corresponds
            !to the value in xc (which is the one corresponding to max loss due to violating feasibility). This leads
            !to crashes/bugs
            return
        end if

        !Compute the loss (unadjusted!!!)
        call check_constr(pars_pntr,c,l,a_pr,rho_pr,aex_pr,constr_ok,loss,.false.,.false.)

        !Now that we know consumption and labour plan, we can compute the current (expected) return
        fc = SWF(c,l,g_st_ind,pars_pntr%P,pars_pntr%mass,pars_pntr%alpha,pars_pntr%A_par,&
        pars_pntr%B_par,pars_pntr%sigma_par,pars_pntr%Gamma_par)

        !Now we have current return in fc, and loss in loss. We just need to compute the continuation return.

        !If the choice is not strictly feasible, force it to be within bounds (extrapolation
        !is very bad, usually leads to divergence). So effectively if the algorithm chooses a point
        !outside of the grid boundaries, we use neares neighbour extapolation (sort of) and add
        !a loss. This is much more stable (otherwise in some corners of the grid the loss if we
        !use extrapolation is compounded and the program diverges). This puts a lower bound on it and it is
        !fine as long as the region where this happens was not optimally reached in the first place (need
        !to experiment with loss).
        if(.not. constr_ok) then
            a_pr=max(pars_pntr%a_min_mat,a_pr)
            a_pr=min(pars_pntr%a_max_mat,a_pr)
            rho_pr=max(pars_pntr%rho_min_vec,rho_pr)
            rho_pr=min(pars_pntr%rho_max_vec,rho_pr)
            aex_pr=max(pars_pntr%aex_min_vec,aex_pr)
            aex_pr=min(pars_pntr%aex_max_vec,aex_pr)
        end if


        !Once we have the point forced within the grid boundaries, get the continuation
        !value by interpolation. We need to get a separate continuation value for all
        !states.


        !Which interpolation function we use depends on the value of VFI_interp_mode (set in parameter file).
        select case(pars_pntr%VFI_interpolation_mode)
            case(1) !Quadrilinear interpolation
                do g_ind = 1,M_par
                    call interp_V(V_old_pntr,[a_pr(g_ind,1),a_pr(g_ind,2),rho_pr(g_ind,1),aex_pr(g_ind,1)],j_guess(g_ind,:)&
                    ,min(g_ind,pars_pntr%M_grid),conval_all(g_ind,1),a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr,size(a1_gr_pntr),&
                    size(rho_gr_pntr),size(aex_gr_pntr))
                end do
            case(2) !Local Shepard's interpolation
                do g_ind = 1,M_par
                    call interp_V_shep(V_old_pntr,[a_pr(g_ind,1),a_pr(g_ind,2),rho_pr(g_ind,1),aex_pr(g_ind,1)],j_guess(g_ind,:)&
                    ,min(g_ind,pars_pntr%M_grid),conval_all(g_ind,1),a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr,size(a1_gr_pntr),&
                    size(rho_gr_pntr),size(aex_gr_pntr),pars_pntr%shep_n,pars_pntr%shep_norm,pars_pntr%shep_norm_val,pars_pntr%shep_floor)
                end do
            case(3) !bspline
                do g_ind = 1,M_par
                    call interp_V_spline(V_old_pntr,[a_pr(g_ind,1),a_pr(g_ind,2),rho_pr(g_ind,1),aex_pr(g_ind,1)],&
                    min(g_ind,pars_pntr%M_grid),conval_all(g_ind,1),a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr,size(a1_gr_pntr),&
                    size(rho_gr_pntr),size(aex_gr_pntr),pars_pntr%bspline_k_a1,pars_pntr%bspline_k_a2,pars_pntr%bspline_k_rho,pars_pntr%bspline_k_aex)
                end do
        end select

        !Now we have continuation value for every possible shock realization. Add the discounted
        !expectation of this term to the current return. (g_st_ind is the last shock realization)

        fc = fc + pars_pntr%beta*maxval(matmul(pars_pntr%P(g_st_ind,:),conval_all))
        !maxval is there simply to convert (1,1) array to scalar.

        !Now fc continues (expected) current plus continuation return for the choice xc

        !Subtract loss stemming from violating some constraints, and multiply by -1 (because we are using a minimization
        !subroutine)
        fc = (fc - loss)*(-1.0_dp)

end subroutine CPC_ret2


!Subroutine CPC_ret_PSO is a copy of CPC_ret but with interface for use with
!particle swarm optimization (NAG FL 24 subroutine E05SAF)

!are passed to this function using common block.
subroutine CPC_ret_PSO(mode,n,xc,fc,vecout,nstate,iuser,ruser)
        integer, intent(inout) :: mode

        integer, intent(in) :: n, nstate
        real (dp), intent(inout) :: fc

        real (dp), intent(inout) :: ruser(*),vecout(n)
        real (dp), intent(in) :: xc(n)
        integer, intent(inout) :: iuser(*)

        real(dp), dimension(:,:,:,:,:), pointer :: V_old_pntr !pointer to V_old, passed using common block
        real(dp), dimension(:), pointer :: a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr !pointers so all grids

        !pointer pars_pntr points to parameters type (it is more convenient to pass it using a common
        !block than create copies of everything we need and put them in a long common block). We
        !just need to be careful to not make any assignment statements because we would then modify
        !parameters elsewhere in the program!!!
        type(par), pointer :: pars_pntr

        !Values of state variables as saved in the grid (hence the suffix st).
        real(dp), dimension(I_par) :: a_st
        real(dp), dimension(I_par-1) :: rho_st
        real(dp) :: aex_st
        integer :: g_st_ind

        !Allocation and next-period states
        real(dp), dimension(M_par,I_par) :: c,l
        real(dp), dimension(M_par,I_par) :: a_pr
        real(dp), dimension(M_par,1) :: rho_pr, aex_pr

        logical :: constr_ok

        integer :: g_ind !index of shock realization

        real(dp) :: loss !loss from violating constraints.

        real(dp), dimension(M_par,1) :: conval_all !continuation value in all states

        integer :: getavfail

        integer, dimension(M_par,4) :: j_guess !The initial guess of position of next-period states
        !in the grid, for evry shock realization.

        !It is very important that nothing that we pass to this function using common block
        !is changed here! Otherwise it might lead to errors elsewhere in the program.
        common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
        !to some subroutines.
        common /pars_pointer/ pars_pntr
        common /grids/ a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr !grids

        common /cpc_ret_bloc/ a_st,rho_st,aex_st,g_st_ind,j_guess

        common /commloss/ loss !For debugging purposes only (so that loss is reported one level higher, not just fc)


        !If mode = 1,2,6,7 - then the program should return first derivatives of objective function.
        !This should never be used (we don't know these and the subroutine should not use local
        !optimizers which require derivatives). If there is an error in this I need to change
        !some optional parameters.
        if(mode == 1 .or. mode == 2 .or. mode == 6 .or. mode == 7) then
        write(*,*) 'ERROR in PSO - wrong value of mode in CPC_ret_PSO. Mode = ',mode
        error stop
        end if

        !For the choice variables xc, get all the remaining variables using subroutine getav
        !Compute all the remaining variables implied by the choice xc
        call get_av(pars_pntr,xc,a_st,rho_st,aex_st,g_st_ind,getavfail,c,l,a_pr,rho_pr,aex_pr,.false.)

        !If non-negativity of consumption or labour is violated, then current return is not defined. Set fc = max_loss
        !(this should guarantee such a point is never deemed optimal) and return.
        if(getavfail > 0) then
            fc = max_loss !maximum loss. This introduces a significant discontinuity, need to be careful about implications elsewhere
            return
        end if

        !Compute the loss (unadjusted!!!)
        call check_constr(pars_pntr,c,l,a_pr,rho_pr,aex_pr,constr_ok,loss,.false.,.false.)

        !Now that we know consumption and labour plan, we can compute the current (expected) return
        fc = SWF(c,l,g_st_ind,pars_pntr%P,pars_pntr%mass,pars_pntr%alpha,pars_pntr%A_par,&
        pars_pntr%B_par,pars_pntr%sigma_par,pars_pntr%Gamma_par)

        !Now we have current return in fc, and loss in loss. We just need to compute the continuation return.

        !If the choice is not strictly feasible, force it to be within bounds (extrapolation
        !is very bad, usually leads to divergence). So effectively if the algorithm chooses a point
        !outside of the grid boundaries, we use neares neighbour extapolation (sort of) and add
        !a loss. This is much more stable (otherwise in some corners of the grid the loss if we
        !use extrapolation is compounded and the program diverges). This puts a lower bound on it and it is
        !fine as long as the region where this happens was not optimally reached in the first place (need
        !to experiment with loss).
        if(.not. constr_ok) then
            a_pr=max(pars_pntr%a_min_mat,a_pr)
            a_pr=min(pars_pntr%a_max_mat,a_pr)
            rho_pr=max(pars_pntr%rho_min_vec,rho_pr)
            rho_pr=min(pars_pntr%rho_max_vec,rho_pr)
            aex_pr=max(pars_pntr%aex_min_vec,aex_pr)
            aex_pr=min(pars_pntr%aex_max_vec,aex_pr)
        end if
        !Once we have the point forced within the grid boundaries, get the continuation
        !value by interpolation. We need to get a separate continuation value for all
        !states.

        !Which interpolation function we use depends on the value of VFI_interp_mode (set in parameter file).
        select case(pars_pntr%VFI_interpolation_mode)
            case(1) !Quadrilinear interpolation
                do g_ind = 1,M_par
                    call interp_V(V_old_pntr,[a_pr(g_ind,1),a_pr(g_ind,2),rho_pr(g_ind,1),aex_pr(g_ind,1)],j_guess(g_ind,:)&
                    ,min(g_ind,pars_pntr%M_grid),conval_all(g_ind,1),a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr,size(a1_gr_pntr),&
                    size(rho_gr_pntr),size(aex_gr_pntr))
                end do
            case(2) !Local Shepard's interpolation
                do g_ind = 1,M_par
                    call interp_V_shep(V_old_pntr,[a_pr(g_ind,1),a_pr(g_ind,2),rho_pr(g_ind,1),aex_pr(g_ind,1)],j_guess(g_ind,:)&
                    ,min(g_ind,pars_pntr%M_grid),conval_all(g_ind,1),a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr,size(a1_gr_pntr),&
                    size(rho_gr_pntr),size(aex_gr_pntr),pars_pntr%shep_n,pars_pntr%shep_norm,pars_pntr%shep_norm_val,pars_pntr%shep_floor)
                end do
            case(3) !bspline
                do g_ind = 1,M_par
                    call interp_V_spline(V_old_pntr,[a_pr(g_ind,1),a_pr(g_ind,2),rho_pr(g_ind,1),aex_pr(g_ind,1)],&
                    min(g_ind,pars_pntr%M_grid),conval_all(g_ind,1),a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr,size(a1_gr_pntr),&
                    size(rho_gr_pntr),size(aex_gr_pntr),pars_pntr%bspline_k_a1,pars_pntr%bspline_k_a2,pars_pntr%bspline_k_rho,pars_pntr%bspline_k_aex)
                end do
        end select

        !Now we have continuation value for every possible shock realization. Add the discounted
        !expectation of this term to the current return. (g_st_ind is the last shock realization)

        fc = fc + pars_pntr%beta*maxval(matmul(pars_pntr%P(g_st_ind,:),conval_all))
        !maxval is there simply to convert (1,1) array to scalar.

        !Now fc continues (expected) current plus continuation return for the choice xc

        !Subtract loss stemming from violating some constraints, and multiply by -1 (because we are using a minimization
        !subroutine)
        fc = (fc - loss)*(-1.0_dp)



end subroutine CPC_ret_PSO


!Subroutine CTS_num_points computes the number of grid points at an iteration of CTS algorithm.
subroutine CTS_num_points(CTS_ind,pars,N_a,N_rho,N_aex)
    integer, intent(in) :: CTS_ind
    type(par), intent(in) :: pars
    integer, intent(out) :: N_a,N_rho,N_aex

    integer :: split = 0 !times that the number of gridpoints is split at this iteration of CTS algorithm

    if(pars%CTS_split_times>0) then
        split = pars%CTS_split_times + 1 - CTS_ind

        N_a = floor(pars%N_a * 0.5_dp**split);
        N_rho = floor(pars%N_rho * 0.5_dp**split);
        N_aex = floor(pars%N_aex * 0.5_dp**split);

        !If some resulting number of grid points is lower than the minimum number of grid points
        !per that variable, adjust this. However, do this only if pars%CTS_split_times > 0. If this is
        !zero it means that we are not using CTS algorithm and it makes no sense to use the floor.
        N_a=maxval([pars%CTS_gridpoints_floor,N_a])
        N_rho=maxval([pars%CTS_gridpoints_floor,N_rho])
        N_aex=maxval([pars%CTS_gridpoints_floor,N_aex])
    else !CTS algorithm not used, use the number of gridpoints given in parameter file
        N_a = pars%N_a
        N_rho = pars%N_rho
        N_aex = pars%N_aex
    end if
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

!subroutine gen_nonequi_grid generates a non-equispaced grid, following
!Heer & Maussner idea.
subroutine gen_nonequi_grid(grid,point1,pointN)
    real(dp), dimension(:) :: grid
    real(dp), intent(in) :: point1,pointN

    real(dp) :: a,b !borders of grid
    integer :: N,N1,N2
    real(dp) :: delta !distance b/w end and final point
    integer :: i

    !If point1 > pointN, it's a mistake (could allow this and reverse the
    !order of points but this will help reveal bugs)
    if(point1>=pointN) then
        write(*,*) 'Error: point1 >= pointN in subroutine gen_equi_grid.'
        error stop
    end if


    N = size(grid) !Total number of points
    N1 = ceiling(real(N)/2.0_dp) !Number of points to be used for the left part of grid
    N2 = (N - N1) + 1 !Number of points to be used for the right part (N1 + N2 = N+1)
    !(so in case of N even, the right part will be slightly more dense)


    !First generate the left part of the grid
    a = point1
    b = (point1 + pointN)/2.0_dp
    delta = (b-a)/(real(N1 - 1,dp)**2.0_dp)
    do i=1,N1
        grid(N1+1-i) = b - ((real(i-1,dp))**2.0_dp)*delta
    end do

    !Now the right part (actually the first point will overwrite an
    !already computed one)
    a = (point1 + pointN)/2.0_dp
    b = pointN
    delta = (b-a)/(real(N2 - 1,dp)**2.0_dp)
    do i=1,N1
        grid(N1+(i-1)) = a + ((real(i-1,dp))**2.0_dp)*delta
    end do

    !do this separately to avoid problems caused by imprecision.
    grid(1) = point1
    grid(N) = pointN

end subroutine gen_nonequi_grid

!Subroutine get_guess_t1 takes as input a solution of the first-period problem for all
!possible shock realizations, and returns a state-contingent plan in the form of the
!choice vector x (as used in simulation of series and VFI).
subroutine gen_guess_t1(FP_result_allm,x_guess_t1)
    type(FP_result_type), dimension(M_par), intent(in) :: FP_result_allm
    real(dp), dimension(ncvar_par) ::  x_guess_t1

    integer :: m !index for cycling over shock realization

    !The first M elements are the labour supply of country 1
    !The following M elements are consumption of country 1
    !The final M-1 elements are the consumption of country 2
    do m = 1,M_par
        x_guess_t1(m) = FP_result_allm(m)%l(1,1)
        x_guess_t1(m+M_par) = FP_result_allm(m)%c(1,1)

        !If M_par > 1, we also save the consumption of country 2,
        !but not for the last shock realization
        if(M_par > 1 .and. m<M_par) then
            x_guess_t1(m+2*M_par) = FP_result_allm(m)%c(1,2)
        end if
    end do

end subroutine gen_guess_t1

!subroutine get_av ('get all variables') computes, for a value of current states, and choice variables x in vector form,
!the value of all relevant variables - the remaining elements of allocation (consumption in country 2 in state M, and labour supply
!in country 2 in all states), and next-period states. It also return the consumption allocation and labour supply in matrix form,
!which is useful in computing return.
!
!If any element of the consumption matrix is <=0, or any element of the labour matrix is < 0, then the current return is undefined
!(for log utility of consumption, otherwise < 0 ). In this case we don't compute the remaining variables, and terminate the
!subroutine with fail = 1 (cons nonpositive), or fail = 2 (labour negative).
subroutine get_av(pars,x,a_st,rho_st,aex_st,g_st_ind,fail,c,l,a_pr,rho_pr,aex_pr,get_cl_only)
    type(par), intent(in) :: pars
    real(dp), dimension(ncvar_par), intent(in) :: x !choice variables in vector form (component of x_guess).

    !Values of current states
    real(dp), dimension(I_par), intent(in) :: a_st
    real(dp), dimension(I_par-1), intent(in) :: rho_st
    real(dp), intent(in) :: aex_st
    integer, intent(in) :: g_st_ind
    logical, intent(in) :: get_cl_only !if true, the subroutine only computes consumption and labour supply.
    !It is useful because sometimes knowing the allocation is all we need (such as when we are checking
    !whether fail > 0, i.e., consumption or labour supply is negative)

    !Output:
    real(dp), dimension(M_par,I_par), intent(out) :: c,l !consumption and labour supply
    integer, intent(out) :: fail !It will be zero if no problems encountered.
    !It will be 1 if consumption was negative (the most serious fail)
    !It will be zero if labour supply was negative (second most serious fail)
    real(dp), dimension(M_par,I_par), intent(out) :: a_pr !next-period adjusted MU-asset holdings
    real(dp), dimension(M_par,1), intent(out) :: rho_pr !next-period ratio of marginal utilities
    real(dp), dimension(M_par,1), intent(out) :: aex_pr !next-period MU-adjusted external debt of country 1

    !Local variables
    real(dp), dimension(M_par,I_par + 1) :: x_mat !choice variables in matrix form
    real(dp), parameter, dimension(M_par,1) :: ones_col = 1.0_dp !column of ones
    real(dp), dimension(M_par,I_par) :: u_c !marginal utility of consumption (element m,i contains MU of consumption of country i in state m)
    real(dp), dimension(1,I_par) :: Eu_c !expected MU of consumption (conditonal on current state) in each country i.

    fail = 0 !initialize

    !The first M elements of x_guess contain labour supply in country 1 in states 1,...,M.
    !Put them into the first column of l
    l(1:M_par,1:1) = reshape(x(1:M_par),[M_par,1])

    !The next M elements contain consumption of country 1 in states 1,...,M
    c(1:M_par,1:1) = reshape(x(1+M_par:2*M_par),[M_par,1])
    !So far we assume I_par = 2 (2 countries). Would need to generalize this for more countries

    !The last M-1 elements contain the consumption of country 2 in states 1,...,M-1
    if(M_par > 1) then
        c(1:M_par-1,2:2) = reshape(x(1+2*M_par:3*M_par-1),[M_par-1,1])
    end if

    !Subroutine c_IM computes the missing element of the consumption matrix, using the constraint equating expected
    !growth of marginal utility
    call c_IM(c,rho_st(1),pars%sigma_par,pars%P,g_st_ind)

    if(c(M_par,I_par) <= 0.0_dp .or. c(M_par,I_par) > 1.0E15_dp) then
        fail = 1 !Some consumption negative or zero
        !if we treat it as zero exactly, then some comparisons elsewhere down the road might not work out).

        !(If utility is logarithmic, we cannot work with 0 consumption, need a strictly positive one).
        !Also - in some perverse case (this requires a very unfortunate combination of transition probabilities,
        !and state variables) the consumption can be infinity, because some terms cancel out and then we take
        !-1 power of 0. While this would result in huge loss from violating other
        !constraints, it causes problems in the NAG optimization subroutines, and
        !it is better to just treat this case separately.
        return
    end if

    !Now the we have consumption in all countries in all states, and labour supply in country 1 in all states, we can compute
    !labour supply in country 2 in all states. It is possible that we get a negative number here (if the aggregate consumption is low,
    !and the production in country 1 is high enough to cover this, plus the government expenditures).

    !Use vectorization to avoid cycles for efficiency
    l(:,2:2) = (pars%mass(1,1)*(c(:,1:1) + pars%G(:,1:1) - pars%theta_0(1,1)*l(:,1:1)) &
                + pars%mass(1,2)*(c(:,2:2) + pars%G(:,2:2)))/(pars%mass(1,2)*pars%theta_0(1,2))

    !Check whether labour supply is less than 0. We only check the labour supply in country 2
    !(this assumes that the input is correct, but that is not always the case. In particular, in PSO,
    !some particles might venture into a territory where this is not true. So rather check all labour supplies.
    !if(any(l(:,2:2)<0.0_dp)) then !(old version)
    if(any(l<0.0_dp)) then
        fail = 2 !Some labour supply smaller than 0 (labour supply equal to zero would not be such an issue)
        return
    end if

    !if get_cl_only, we only care about consumption and labour. Set everything else to zero
    !to avoid potential error, and then return.
    if(get_cl_only) then
        a_pr = 0.0_dp
        rho_pr = 0.0_dp
        aex_pr = 0.0_dp
        return
    end if

    !Now that we have consumption and labour supply, compute next-period states.

    !The expected marginal utility of consumption is used repeatedly. In the first step compute that for both countries and
    !save it (row vector)
    !Get marginal utilities of all countries in all states, save these in matrix u_c
    call util_c(c,u_c,pars%A_par,pars%sigma_par)
    !Expected utility (1 x I row vector)
    Eu_c = matmul(pars%P(g_st_ind:g_st_ind,:),u_c)

    !Compute next-period ratio of marginal utilities for all states. Return this as Mx1 column vector
    rho_pr(:,1) = u_c(:,1)/u_c(:,2)

    !Compute next-period assets (total debt of both countries). These will depend on functional form of utility function. Defined in module mod_util
    !a_pr,aex_pr
    call a_prime(a_st,c,l,Eu_c,pars%mass,pars%theta_0,pars%G,pars%B_par,pars%sigma_par,pars%Gamma_par,pars%beta,a_pr)
    call aex_prime(aex_st,c,l,Eu_c,pars%mass,pars%theta_0,pars%G,pars%sigma_par,pars%beta,aex_pr)

end subroutine get_av

subroutine get_av_debug(pars,x,a_st,rho_st,aex_st,g_st_ind,fail,c,l,a_pr,rho_pr,aex_pr,get_cl_only)
    type(par), intent(in) :: pars
    real(dp), dimension(ncvar_par), intent(in) :: x !choice variables in vector form (component of x_guess).

    !Values of current states
    real(dp), dimension(I_par), intent(in) :: a_st
    real(dp), dimension(I_par-1), intent(in) :: rho_st
    real(dp), intent(in) :: aex_st
    integer, intent(in) :: g_st_ind
    logical, intent(in) :: get_cl_only !if true, the subroutine only computes consumption and labour supply.
    !It is useful because sometimes knowing the allocation is all we need (such as when we are checking
    !whether fail > 0, i.e., consumption or labour supply is negative)

    !Output:
    real(dp), dimension(M_par,I_par), intent(out) :: c,l !consumption and labour supply
    integer, intent(out) :: fail !It will be zero if no problems encountered.
    !It will be 1 if consumption was negative (the most serious fail)
    !It will be zero if labour supply was negative (second most serious fail)
    real(dp), dimension(M_par,I_par), intent(out) :: a_pr !next-period adjusted MU-asset holdings
    real(dp), dimension(M_par,1), intent(out) :: rho_pr !next-period ratio of marginal utilities
    real(dp), dimension(M_par,1), intent(out) :: aex_pr !next-period MU-adjusted external debt of country 1

    !Local variables
    real(dp), dimension(M_par,I_par + 1) :: x_mat !choice variables in matrix form
    real(dp), parameter, dimension(M_par,1) :: ones_col = 1.0_dp !column of ones
    real(dp), dimension(M_par,I_par) :: u_c !marginal utility of consumption (element m,i contains MU of consumption of country i in state m)
    real(dp), dimension(1,I_par) :: Eu_c !expected MU of consumption (conditonal on current state) in each country i.

    fail = 0 !initialize

    !The first M elements of x_guess contain labour supply in country 1 in states 1,...,M.
    !Put them into the first column of l
    l(1:M_par,1:1) = reshape(x(1:M_par),[M_par,1])

    !The next M elements contain consumption of country 1 in states 1,...,M
    c(1:M_par,1:1) = reshape(x(1+M_par:2*M_par),[M_par,1])
    !So far we assume I_par = 2 (2 countries). Would need to generalize this for more countries

    !The last M-1 elements contain the consumption of country 2 in states 1,...,M-1
    if(M_par > 1) then
        c(1:M_par-1,2:2) = reshape(x(1+2*M_par:3*M_par-1),[M_par-1,1])
    end if

    write(*,*) 'g_St_ind in get_av_debug b4 calling c_IM_debug = ',g_st_ind

    !Subroutine c_IM computes the missing element of the consumption matrix, using the constraint equating expected
    !growth of marginal utility
    call c_IM_debug(c,rho_st(1),pars%sigma_par,pars%P,g_st_ind)

    if(c(M_par,I_par) <= 0.0_dp .or. c(M_par,I_par) > 1.0E15) then
        fail = 1 !Some consumption negative or zero
        return
    end if



    !Now the we have consumption in all countries in all states, and labour supply in country 1 in all states, we can compute
    !labour supply in country 2 in all states. It is possible that we get a negative number here (if the aggregate consumption is low,
    !and the production in country 1 is high enough to cover this, plus the government expenditures).

    !Use vectorization to avoid cycles for efficiency
    l(:,2:2) = (pars%mass(1,1)*(c(:,1:1) + pars%G(:,1:1) - pars%theta_0(1,1)*l(:,1:1)) &
                + pars%mass(1,2)*(c(:,2:2) + pars%G(:,2:2)))/(pars%mass(1,2)*pars%theta_0(1,2))

    !Check whether labour supply is less than 0. We only check the labour supply in country 2
    if(any(l(:,2:2)<0.0_dp)) then
        fail = 2 !Some labour supply smaller than 0 (labour supply equal to zero would not be such an issue)
        return
    end if

    !if get_cl_only, we only care about consumption and labour. Set everything else to zero
    !to avoid potential error, and then return.
    if(get_cl_only) then
        a_pr = 0.0_dp
        rho_pr = 0.0_dp
        aex_pr = 0.0_dp
        return
    end if

    !Now that we have consumption and labour supply, compute next-period states.

    !The expected marginal utility of consumption is used repeatedly. In the first step compute that for both countries and
    !save it (row vector)
    !Get marginal utilities of all countries in all states, save these in matrix u_c
    call util_c(c,u_c,pars%A_par,pars%sigma_par)
    !Expected utility (1 x I row vector)
    Eu_c = matmul(pars%P(g_st_ind:g_st_ind,:),u_c)


    write(*,*) 'debug in get_av_debug:'
    write(*,*) 'c = ',c
    write(*,*) 'l = ',l
    write(*,*) 'u_c = ',u_c

    write(*,*) 'rhoprman_st1 = ',u_c(1,1)/u_c(1,2)
    write(*,*) 'rhoprman_st2 = ',u_c(2,1)/u_c(2,2)




    !Compute next-period ratio of marginal utilities for all states. Return this as Mx1 column vector
    rho_pr(:,1) = u_c(:,1)/u_c(:,2)

    !rho_pr is displayed outside of this subroutine anyway
    write(*,*) '_______________________________________'

    !Compute next-period assets (total debt of both countries). These will depend on functional form of utility function. Defined in module mod_util
    !a_pr,aex_pr
    call a_prime(a_st,c,l,Eu_c,pars%mass,pars%theta_0,pars%G,pars%B_par,pars%sigma_par,pars%Gamma_par,pars%beta,a_pr)
    call aex_prime(aex_st,c,l,Eu_c,pars%mass,pars%theta_0,pars%G,pars%sigma_par,pars%beta,aex_pr)

end subroutine get_av_debug


!subroutine get_av_fp is pretty much a copy of get_av. The difference is that it gives us the values
!of all variables for choice variables in the first period, when we do not have the rho state (and
!constraint). Also, we only choose variables in one state. The choice is consumption of both countries,
!and labor supply in country 1. There are therefore 3 choice variables (needs to be generalized
!for I > 2). Vector x contains the choice. x(1) is labour supply of country 1, x(2) is consumption in
!country 1, x(3) is consumption in country 2. Also, the states are initial conditions.
subroutine get_av_fp(pars,x,b_init,b_ex_init,g_st_ind,fail,c,l,b_pr,b_ex_pr,a_pr,rho_pr,aex_pr,get_cl_only)
    type(par), intent(in) :: pars
    real(dp), dimension(3), intent(in) :: x !choice variables in vector form (component of x_guess).

    !Values of current states
    real(dp), dimension(I_par), intent(in) :: b_init
    real(dp), intent(in) :: b_ex_init
    integer, intent(in) :: g_st_ind !index of initial state
    logical, intent(in) :: get_cl_only !if true, the subroutine only computes consumption and labour supply.
    !It is useful because sometimes knowing the allocation is all we need (such as when we are checking
    !whether fail > 0, i.e., consumption or labour supply is negative)

    !Output:
    real(dp), dimension(1,I_par), intent(out) :: c,l !consumption and labour supply
    integer, intent(out) :: fail !It will be zero if no problems encountered.
    !It will be 1 if consumption was negative (the most serious fail)
    !It will be zero if labour supply was negative (second most serious fail)
    real(dp), dimension(1,I_par), intent(out) :: a_pr!next-period adjusted MU-asset holdings
    real(dp), dimension(1,1), intent(out) :: rho_pr !next-period ratio of marginal utilities
    real(dp), dimension(1,1), intent(out) :: aex_pr !next-period MU-adjusted external debt of country 1

    !Also return b_pr and b_ex_pr. Sometimes it is useful to get this, and we compute this anyway.
    real(dp), dimension(1,1), intent(out) :: b_ex_pr
    real(dp), dimension(1,I_par), intent(out) :: b_pr


    real(dp), dimension(1,I_par) :: u_c !marginal utility of consumption (element m,i contains MU of consumption of country i in state m)

    fail = 0 !initialize
    l = 0.0_dp

    !First get the input into elements of the consumption and labour vector
    l(1,1) = x(1) !labour supply of country 1
    c(1,1) = x(2) !consumption of country 1
    c(1,2) = x(3) !labour supply of country 1

    !Get consumption of country 2.
    l(:,2:2) = (pars%mass(1,1)*(c(:,1:1) + pars%G(g_st_ind:g_st_ind,1:1) - pars%theta_0(1,1)*l(:,1:1)) &
                + pars%mass(1,2)*(c(:,2:2) + pars%G(g_st_ind:g_st_ind,2:2)))/(pars%mass(1,2)*pars%theta_0(1,2))

    !Check whether labour or consumption supply is less than 0. Technically it should never happen that we pass
    !some x containing a negative element, so we wouldn't need to check all but only labour supply in country 2,
    !but do this just in case. It is not costly and could prevent serious bugs.

    if(any(c<0.0_dp)) then
        fail = 1 !Some labour supply smaller than 0 (labour supply equal to zero would not be such an issue)
        return
    end if

    if(any(l<0.0_dp)) then
        fail = 2 !Some labour supply smaller than 0 (labour supply equal to zero would not be such an issue)
        return
    end if


    !if get_cl_only, we only care about consumption and labour. Set everything else to zero
    !to avoid potential error, and then return.
    if(get_cl_only) then
        a_pr = 0.0_dp
        rho_pr = 0.0_dp
        aex_pr = 0.0_dp
        return
    end if


    !Now that we have consumption and labour supply, compute next-period states.
    !Get marginal utilities of all countries in all states, save these in matrix u_c
    call util_c(c,u_c,pars%A_par,pars%sigma_par)

    !Compute next-period ratio of marginal utilities for all states. Return this as Mx1 column vector
    rho_pr(:,1) = u_c(:,1)/u_c(:,2)

    !Compute next-period assets (total debt of both countries). These depend on functional form so for consistency,
    !are defined in subroutine b_prime.
    call b_prime_fp(b_init,c,l,pars%G(g_st_ind:g_st_ind,1:I_par),pars%mass,pars%theta_0,pars%B_par,pars%sigma_par,pars%gamma_par,b_pr)
    call b_ex_prime_fp(b_ex_init,c,l,pars%G(g_st_ind:g_st_ind,1:I_par),pars%mass,pars%theta_0,b_ex_pr)

    !Finally, multiply these by marginal utility and divide by beta to get states a_pr,aex_pr
    a_pr = b_pr*u_c/pars%beta
    aex_pr = b_ex_pr*u_c(1,1)/pars%beta

end subroutine get_av_fp

!subroutine get_aut_ss computes a deterministic steady state in in autarky (isolated economies)
!with asset holdings equal to zero. The government expenditure in a country is the expected
!government expenditure (according to steady-state probabilities) maintained forever if s_ind = -1.
!Otherwise the government expenditure is forever maintained at the level corresponding to s_ind.
!The latter is used for getting initial guess in the first period (it is better because of symmetry,
!at least if asset holdings are initially zero).
subroutine get_aut_ss(pars,SS,s_ind)
    type(par), intent(in) :: pars
    type(SS_type), intent(out) :: SS
    integer, intent(in) :: s_ind

    real(dp) :: tau_ss(I_par) !temporary ss tax rate


    !For now, this is implemented only for log utility of consumption (sigma = 1). Otherwise we need to use an iterative procedure which is fine,
    !But we need to make sure that we rule out the solution in which consumption is zero.
    if(pars%sigma_par /= 1) then
        !Try to get the steady state using the previously written subroutine get_ss. This struggles
        !with non-unique solution and sometimes finding a zero consumption solution,
        !usually when initial assets are zero. This one actually takes into account the initial assets. This is
        !perhaps a mistake and it should eb

        !The steady state now actually is not an autarky steady state - but as long as initial assets are zero it is.
        !the same. I should fix this later.
        call get_ss(pars,SS,s_ind)

        !If it fails (we get zeros) then try to use the logarithmic utility steady state (but that might be
        !quite bad).


        if(this_image() == 1) then
        write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        write(*,*) 'Warning: Using subroutine get_ss to get autarky steady state. This might be problematic.'
        !Hopefully it will not be but write out the consumption just in case and also compute the steady state tax
        !rate in both countries and display it.
        end if


    if(this_image() == 1) then
        write(*,*) 'c_ss = ',SS%c
        write(*,*) 'l_ss = ',SS%l
        !compute steady state tax rate
        if(I_par == 2) then
            tau_ss(1) = 1.0_dp - (pars%B_par/pars%theta_0(1,1))*(ss%c(1,1)**pars%sigma_par)*(ss%l(1,1)**pars%gamma_par)
            tau_ss(2) = 1.0_dp - (pars%B_par/pars%theta_0(1,2))*(ss%c(1,2)**pars%sigma_par)*(ss%l(1,2)**pars%gamma_par)
            write(*,*) 'tau_ss = ',tau_ss
        end if
    end if

    !If all the elements are significantly far from zero, then return this. Otherwise do not
    !in which case the logarithmic SS below will be used.
    if(SS%c(1,1) > 0.001_dp .and. SS%c(1,2) > 0.001_dp .and.SS%l(1,1) > 0.001_dp .and.SS%l(1,2) > 0.001_dp) then
        return
    else
        if(this_image() ==1) write(*,*) 'Some elements are too close to zero, using logarithmic steady state instead'
    end if

    end if


    !If we got to this part of the program, we can assume logarithmic utility of consumption

    !IF the input parameter s_ind = -1, use the expected government expenditure. Otherwise use the government
    !expenditure which corresponds to s_ind
    if(s_ind == -1) then
        SS%l = (1.0_dp/pars%B_par)**(1.0_dp/(1.0_dp+pars%gamma_par))
        SS%c = SS%l* pars%theta_0 - reshape(pars%Eg,[1,I_par])
        !total consumption
        SS%c_ss_tot = pars%mass(1,1)*SS%c(1,1) + pars%mass(1,2)*SS%c(1,2)
        !Excess of government expenditures in country 2 over country 1, scaled by LFFC_kappa
        SS%kappa_g2ming1 = (pars%Eg(1,2) - pars%Eg(1,1))*pars%LFFC_kappa/2.0_dp
        !divide by 2 so that if kappa = 1, and we subtract and add this to consumption,
        !we get a difference in consumption which corresponds to the difference in
        !government expenditure - not double of that. IT doesn't really matter, it's just
        !so that kappa = 1 is the benchmark of no smoothing.
    else
        SS%l = (1.0_dp/pars%B_par)**(1.0_dp/(1.0_dp+pars%gamma_par))
        SS%c = SS%l* pars%theta_0 - pars%G(s_ind:s_ind,:)
        !total consumption
        SS%c_ss_tot = pars%mass(1,1)*SS%c(1,1) + pars%mass(1,2)*SS%c(1,2)
        !Excess of government expenditures in country 2 over country 1, scaled by LFFC_kappa
        SS%kappa_g2ming1 = (pars%G(s_ind,2) - pars%G(s_ind,1))*pars%LFFC_kappa/2.0_dp
    end if

    !Also in the case of logarithmic utility display information about steady state
    if(this_image() == 1) then
        write(*,*) 'c_ss = ',SS%c
        write(*,*) 'l_ss = ',SS%l
        !compute steady state tax rate
        if(I_par == 2) then
            tau_ss(1) = 1.0_dp - (pars%B_par/pars%theta_0(1,1))*(ss%c(1,1)**pars%sigma_par)*(ss%l(1,1)**pars%gamma_par)
            tau_ss(2) = 1.0_dp - (pars%B_par/pars%theta_0(1,2))*(ss%c(1,2)**pars%sigma_par)*(ss%l(1,2)**pars%gamma_par)
            write(*,*) 'tau_ss = ',tau_ss
        end if
    end if

    !Otherwise (more general preferences) we would have to solve a non-linear equation of 1 variable (and rule out the zero consumption solution)
end subroutine get_aut_ss

!subroutine get_guess_img computes an initial guess
!at every gridpoint on which an image works. It also computes the share of nodes at which this is feasible.
!(so that the value of loss function is zero).
!This result is saved in coarray share_feas_img.
!The feasible guesses are returned in vector of type c_guess_all (see its definition in this module).
subroutine get_guess_img(x_guess_img,SS,min_ind_img,max_ind_img,ind_unfold_all,pars,grids,share_feas_img,&
avg_loss_img,CTS_ind,x_pol_all_unf,grids_tmp)
    real(dp), codimension[*] :: share_feas_img !The weighted share of points where a feasible c was found
    real(dp) :: avg_loss_img
    !type(SS_type), intent(in) :: SS

    type(ss_type), dimension(0:M_par), intent(in) :: SS !steady-state. The first element (index 0) corresponds to steady state when
    !government expenditure is kept forever at the expected (long-run SS probabilities) government expenditure.
    !Element with index m contains the same, but the expenditure corresponds to shock realization with index m.
    !This is useful for getting initial guess in first period, and for debugging.

    integer, intent(in) :: min_ind_img,max_ind_img
    type(par), codimension[*] :: pars !parameters. Not sure why this is a coarray here... Could be a coarray
    !only where we need to share data, here it could be just a standard type. Low priority issue though.
    real(dp), dimension(ncvar_par,min_ind_img:max_ind_img), intent(out) :: x_guess_img !x_guess_img
    !superceded c_guess_img (it contains choice of both consumption and labour supply, See defin. of type x_guess).
    integer, dimension(5,min_ind_img:max_ind_img) :: ind_unfold_all
    type(grids_type), intent(in) :: grids
    integer, intent(in) :: CTS_ind
    real(dp), dimension(:,:,:,:,:,:), allocatable, intent(in) :: x_pol_all_unf !Policy function from the previous stage of CTS
    type(grids_type), intent(in) :: grids_tmp

    integer :: maxatt_iter

    integer :: gp_ind !index for cycling over gridpoints. The main index from which we will get the indices
    !for the state variables from ind_unfold_all
    !Some of the following variables will become row vectors.
    !integer :: a1_ind,a2_ind,rho_ind,theta_ind,m_ind !- indices not needed to be kept separately?

    !Values of state variables as saved in the grid (hence the suffix st).
    real(dp), dimension(I_par) :: a_st
    real(dp), dimension(I_par-1) :: rho_st
    real(dp) :: aex_st,g_st
    integer :: g_st_ind

    real(dp) :: loss_gp, loss_img !loss at gridpoint, loss at image (average weighted by number of
    !gridpoints at this image).

    if(I_par>2) then
        write(*,*) 'Error: subroutine get_guess_img assumes I_par = 2.'
        error stop
    end if


    !If this is the first iteration of CTS algorithm (and it is not the only iteration) then
    !the number of grid points is fairly low. We multiply the number of attempts in this iteration
    !of CTS tenfold (somewhat arbitrary). The number of attempts actually makes a big difference.
    if(CTS_ind == 1 .and. pars%CTS_split_times > 0) then
        maxatt_iter = pars%LFFC_attmax * 10
    else
        maxatt_iter = pars%LFFC_attmax
    end if


    !Generate an initial guess using the steady-state
    !and save it in all elements of x_guess_img. These will then be passed to subroutine LFFC (intent(inout))
    !and overwritten on output. First get this guess once and save into the element of x_guess_img associated with
    !the first gridpoint that this image works on, then copy it to all other elements.
    ![Note: this assumes I_par = 2]

    !Actually, the way the program currently works is that LFFC only uses the steady-state consumption
    !of both countries in all states to compute the total consumption in a state (this is constant if
    !M_par is not equal to 2, in which case we use an average level of consumption).


    !In the new version of the program, I use (for each SS), a steady-state in autarky
    !which would be a SS if the government expenditure was kept at that level forever.
    !This improves the guess mainly for the case of agregate shocks.
    !The adjustment only works for M_par == 2, otherwise an old version of the program is used
    if(M_par == 2) then
    !Fill in the 5 elements by hand. Firs of use SS(0) to make sure that nothing is changed in this
    !case.
        !The first M elements are the labour supply of country 1 (in autarky ss)
        x_guess_img(1,min_ind_img:max_ind_img) = SS(0)%l(1,1)
        x_guess_img(2,min_ind_img:max_ind_img) = SS(0)%l(1,1)
        !The next M elements are steady state consumption of country 1
        x_guess_img(3,min_ind_img:max_ind_img) = SS(0)%c(1,1)
        x_guess_img(4,min_ind_img:max_ind_img) = SS(0)%c(1,1)
        !The final M-1 elements are the SS consumption of country 2
        x_guess_img(5,min_ind_img:max_ind_img) = SS(0)%c(1,2)
    else
        !The first M elements are the labour supply of country 1 (in autarky ss)
        x_guess_img(1:M_par,min_ind_img:max_ind_img) = SS(0)%l(1,1)
        !The next M elements are steady state consumption of country 1
        x_guess_img(1+M_par:2*M_par,min_ind_img:max_ind_img) = SS(0)%c(1,1)
        !The final M-1 elements are the SS consumption of country 2
        if(M_par > 1) then
            x_guess_img(1+2*M_par:3*M_par-1,min_ind_img:max_ind_img) = SS(0)%c(1,2)
        end if
    end if


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
        aex_st = grids%aex_gr(ind_unfold_all(4,gp_ind))
        g_st_ind = ind_unfold_all(5,gp_ind)

        call LFFC(a_st,rho_st,aex_st,g_st_ind,pars,x_guess_img(1:ncvar_par,gp_ind),loss_gp,CTS_ind,x_pol_all_unf,grids_tmp,SS,maxatt_iter)

        !update statistics on feasibility of shares, etc.
        avg_loss_img = avg_loss_img + loss_gp/real(max_ind_img-min_ind_img+1)

        !If loss is less than 0.001, then treat the point as strictly feasible for debugging purposes.
        if(loss_gp < 0.001_dp) then
            share_feas_img = share_feas_img + 1.0_dp/real(max_ind_img-min_ind_img+1)
        end if

    end do
end subroutine get_guess_img


!subroutine get_ss computes a deterministic steady state in which asset holdings are
!kept constant at the initial level. The government expenditure in a country is the expected
!government expenditure (according to steady-state probabilities) maintained forever.
!The issue with this is that the steady-state seems to be non-unique and quite often
!we get a steady state with 0 consumption, which is terrible as an initial guess.
!Therefore we use the autarky steady state instead. This subroutine is left here
!only as a work in progress (it could be useful in the future).
subroutine get_ss(pars,SS,s_ind)
    type(par), intent(in) :: pars
    type(SS_type), intent(out) :: SS
    integer :: s_ind

    type(par) :: pars_copy
    real(dp), dimension(1,I_par) :: Eg !expected government expenditures in all countries

    !variables used by NAG to solve the system of equations
    integer, parameter :: n = 4 !number of equations
    real(dp), dimension(n) :: x !choice variables vector (c1,c2,n1,n2)
    real(dp), dimension(n) :: fvec !values of the n equations for given x
    real(dp), dimension(n,n) :: fjac !jacobian
    real(dp) :: x_tol = 0.00000000001_dp
    integer, dimension(1) :: iuser = 0
    real(dp), dimension(1) :: ruser = 0
    integer :: ifail

    integer :: s_ind_copy

    common pars_copy, s_ind_copy !a copy of parameters which will be put in a common block so that it can be accessed
    !by subroutine ss_sys (steady-state system)

    !create copies for common block
    pars_copy = pars
    s_ind_copy = s_ind



    !Use the solver (NAG subroutine C05RBF)
    ifail = -1 !no stopping but display warning message in case of fail
    x = 1.0_dp !start with all consumptions and labour supplies equal to 1.

!    !Using Jacobian - this one still was not fixed (a few typos in signs - but the one without Jacobian works just fine).
!    call C05RBF(ss_sys, n, x, fvec, fjac, x_tol, iuser, ruser, ifail)

    !Not using Jacobian
    call C05QBF(ss_sys2, n, x, fvec, x_tol, iuser, ruser, ifail)


    !Now that we have the solution, copy it to SS_type
    SS%c(1,1) = x(1)
    SS%c(1,2) = x(2)
    SS%l(1,1) = x(3)
    SS%l(1,2) = x(4)


end subroutine get_ss

!Subroutine ss_sys is the steady state system of equations in form usable by NAG
subroutine ss_sys(N, X, FVEC, FJAC, IUSER, RUSER, IFLAG)
    integer :: N, IUSER(*), IFLAG
    real(dp) :: X(N), FVEC(N), FJAC(N,N), RUSER(*)

    integer :: s_ind_copy
    type(par) :: pars_copy

    common pars_copy, s_ind_copy

    !as government expenditure for steady state computation, use the expected government expenditure
    !in both countries (using steady-state probabilities). This is computed in subroutine read_parameters
    !in mod_pars (see there for details) and saved in pars_copy%Eg (1,I_par vector)

    !Remember: c^1 = x(1), c^2 = x(2), n^1 = x(3), n^2 = x(4) [x=(c^1,c^2,n^1,n^2)]

    !The following still contains an error in which there was a wrong sign.
    !Basically it should be +pars_copy%B_par* (x(1)**pars_copy%sigma_par), not minus. And elsewhere this
    !error is repeated.

    !Compute the function value (in f(x) = 0) equation
    fvec(1) = pars_copy%mass(1,1) * (pars_copy%Eg(1,1) - pars_copy%theta_0(1,1)*x(3) - &
        pars_copy%B_par* (x(1)**pars_copy%sigma_par) * (x(3)**(pars_copy%gamma_par + 1.0_dp)) ) - &
        ((pars_copy%beta - 1.0_dp)/pars_copy%beta)*pars_copy%b_init(1,1)

    fvec(2) = pars_copy%mass(1,2) * (pars_copy%Eg(1,2) - pars_copy%theta_0(1,2)*x(4) - &
        pars_copy%B_par* (x(2)**pars_copy%sigma_par) * (x(4)**(pars_copy%gamma_par + 1.0_dp)) ) - &
        ((pars_copy%beta - 1.0_dp)/pars_copy%beta)*pars_copy%b_init(1,2)

    fvec(3) = pars_copy%mass(1,1) * (x(1) + pars_copy%Eg(1,1) - pars_copy%theta_0(1,1)*x(3)) + &
        pars_copy%mass(1,2) * (x(2) + pars_copy%Eg(1,2) - pars_copy%theta_0(1,2)*x(4))

    fvec(4) = pars_copy%mass(1,1) * (x(1) + pars_copy%Eg(1,1) - pars_copy%theta_0(1,1)*x(3)) - &
        ((pars_copy%beta - 1.0_dp)/pars_copy%beta)*pars_copy%b_ex_init


    !Compute the Jacobian matrix (only if IFLAG == 2 - see documentation of the NAG subroutine)
    if(iflag == 2) then
        fjac(1,1) = -pars_copy%mass(1,1)*pars_copy%B_par*pars_copy%sigma_par* (x(1)** (pars_copy%sigma_par - 1.0_dp)) * &
            (x(3)**(1.0_dp + pars_copy%gamma_par))
        fjac(1,2) = 0.0_dp
        fjac(1,3) = -pars_copy%mass(1,1)*(pars_copy%theta_0(1,1) + pars_copy%B_par*(1.0_dp + pars_copy%gamma_par) * &
            (x(1)**(pars_copy%sigma_par)) * (x(3)**(pars_copy%gamma_par)) )
        fjac(1,4) = 0.0_dp

        fjac(2,1) = 0.0_dp
        fjac(2,2) = -pars_copy%mass(1,2)*pars_copy%B_par*pars_copy%sigma_par* (x(2)** (pars_copy%sigma_par - 1.0_dp)) * &
            (x(4)**(1.0_dp + pars_copy%gamma_par))
        fjac(2,3) = 0.0_dp
        fjac(2,4) = -pars_copy%mass(1,2)*(pars_copy%theta_0(1,2) + pars_copy%B_par*(1.0_dp + pars_copy%gamma_par) * &
            (x(2)**(pars_copy%sigma_par)) * (x(4)**(pars_copy%gamma_par)) )

        fjac(3,1) = pars_copy%mass(1,1)
        fjac(3,2) = pars_copy%mass(1,2)
        fjac(3,3) = -pars_copy%mass(1,1)*pars_copy%theta_0(1,1)
        fjac(3,4) = -pars_copy%mass(1,2)*pars_copy%theta_0(1,2)

        fjac(4,1) = pars_copy%mass(1,1)
        fjac(4,2) = 0.0_dp
        fjac(4,3) = -pars_copy%mass(1,1)*pars_copy%theta_0(1,1)
        fjac(4,4) = 0.0_dp
    end if

end subroutine ss_sys

!Subroutine ss_sys is a copy of the above with a different interface so it can be used by a different
!optimization subroutine.
subroutine ss_sys2(N, X, FVEC, IUSER, RUSER, IFLAG)
    integer :: N, IUSER(*), IFLAG
    real(dp) :: X(N), FVEC(N), RUSER(*)

    integer :: s_ind_copy
    type(par) :: pars_copy

    common pars_copy,s_ind_copy

    !as government expenditure for steady state computation, use the expected government expenditure
    !in both countries (using steady-state probabilities). This is computed in subroutine read_parameters
    !in mod_pars (see there for details) and saved in pars_copy%Eg (1,I_par vector)

    !Remember: c^1 = x(1), c^2 = x(2), n^1 = x(3), n^2 = x(4) [x=(c^1,c^2,n^1,n^2)]

    !Compute the function value (in f(x) = 0) equation

    !if(s_ind = -1, use the expected value of government expenditures). Otherwise use the current value
    !corresponding to s_ind

    if(s_ind_copy == -1) then
    fvec(1) = pars_copy%mass(1,1) * (pars_copy%Eg(1,1) - pars_copy%theta_0(1,1)*x(3) + &
        pars_copy%B_par* (x(1)**pars_copy%sigma_par) * (x(3)**(pars_copy%gamma_par + 1.0_dp)) ) - &
        ((pars_copy%beta - 1.0_dp)/pars_copy%beta)*pars_copy%b_init(1,1)

    fvec(2) = pars_copy%mass(1,2) * (pars_copy%Eg(1,2) - pars_copy%theta_0(1,2)*x(4) + &
        pars_copy%B_par* (x(2)**pars_copy%sigma_par) * (x(4)**(pars_copy%gamma_par + 1.0_dp)) ) - &
        ((pars_copy%beta - 1.0_dp)/pars_copy%beta)*pars_copy%b_init(1,2)

    fvec(3) = pars_copy%mass(1,1) * (x(1) + pars_copy%Eg(1,1) - pars_copy%theta_0(1,1)*x(3)) + &
        pars_copy%mass(1,2) * (x(2) + pars_copy%Eg(1,2) - pars_copy%theta_0(1,2)*x(4))

    fvec(4) = pars_copy%mass(1,1) * (x(1) + pars_copy%Eg(1,1) - pars_copy%theta_0(1,1)*x(3)) - &
        ((pars_copy%beta - 1.0_dp)/pars_copy%beta)*pars_copy%b_ex_init
    else

    fvec(1) = pars_copy%mass(1,1) * (pars_copy%g(s_ind_copy,1) - pars_copy%theta_0(1,1)*x(3) + &
        pars_copy%B_par* (x(1)**pars_copy%sigma_par) * (x(3)**(pars_copy%gamma_par + 1.0_dp)) ) - &
        ((pars_copy%beta - 1.0_dp)/pars_copy%beta)*pars_copy%b_init(1,1)

    fvec(2) = pars_copy%mass(1,2) * (pars_copy%g(s_ind_copy,2) - pars_copy%theta_0(1,2)*x(4) + &
        pars_copy%B_par* (x(2)**pars_copy%sigma_par) * (x(4)**(pars_copy%gamma_par + 1.0_dp)) ) - &
        ((pars_copy%beta - 1.0_dp)/pars_copy%beta)*pars_copy%b_init(1,2)

    fvec(3) = pars_copy%mass(1,1) * (x(1) + pars_copy%g(s_ind_copy,1) - pars_copy%theta_0(1,1)*x(3)) + &
        pars_copy%mass(1,2) * (x(2) + pars_copy%g(s_ind_copy,2) - pars_copy%theta_0(1,2)*x(4))

    fvec(4) = pars_copy%mass(1,1) * (x(1) + pars_copy%g(s_ind_copy,1) - pars_copy%theta_0(1,1)*x(3)) - &
        ((pars_copy%beta - 1.0_dp)/pars_copy%beta)*pars_copy%b_ex_init
    end if

end subroutine ss_sys2


!subroutine for generating grids. 4 equispaced grids are generated for the first four states
!(a_1,a_2,rho,t).
!The grid for shock space is loaded from parameters and we keep it separate (this grid is
!never resized).
subroutine grids_gen(A1_gr,A2_gr,rho_gr,aex_gr,N_a,N_rho,N_aex,pars)
    real(dp), dimension(:), intent(out) :: A1_gr,A2_gr,rho_gr,aex_gr
    integer, intent(in) :: N_a,N_rho,N_aex !can't just load these from parameter type
    !because in CTS algorithm they vary.
    type(par) :: pars

    integer :: i,rho_mid_ind !(index for cycling over grid, middle index of the grid for rho)

    !if grids_type is 1 (all grids equispaced) or 3 (all grids equispaced
    !except for grid for rho), generate the equispaced grids:
    if(pars%grids_type == 1 .or. pars%grids_type == 3) then
        call gen_equi_grid(A1_gr,pars%a_min(1,1),pars%a_max(1,1))
        call gen_equi_grid(A2_gr,pars%a_min(1,2),pars%a_max(1,2))
        call gen_equi_grid(aex_gr,pars%aex_min,pars%aex_max)
    else
        call gen_nonequi_grid(A1_gr,pars%a_min(1,1),pars%a_max(1,1))
        call gen_nonequi_grid(A2_gr,pars%a_min(1,2),pars%a_max(1,2))
        call gen_nonequi_grid(aex_gr,pars%aex_min,pars%aex_max)
    end if

    !for rho_grid, only grids_type == 1 corresponds to equispaced grid.
    !In addition, if we use symmetry, we want to make sure that the middle of the grids
    !is 1 (so that 1 is in the grid).
    if(pars%rho_gr_symm == .false.) then
    !in this case proceed as for the other grids
        if(pars%grids_type == 1) then
            call gen_equi_grid(rho_gr,pars%rho_min,pars%rho_max)
        else
            call gen_nonequi_grid(rho_gr,pars%rho_min,pars%rho_max)
        end if
    else
        if(pars%grids_type == 1) then
            call gen_equi_grid(rho_gr,pars%rho_min,2.0_dp-pars%rho_min)
        else
            call gen_nonequi_grid(rho_gr,pars%rho_min,2.0_dp-pars%rho_min)
        end if
        !Now set the 'right half' of the grid (points to the right of 1) using the
        !other half (so that there is a corresponding inverted point for every point on the left of 1).
        !We already know that the number of grid points is odd because that's set earlier in the program
        !whenever rho_gr_symm = .true.
        rho_mid_ind = N_rho/2 + 1
        do i = 1,(N_rho/2)
            rho_gr(rho_mid_ind + i) = 1.0_dp/rho_gr(rho_mid_ind - i)
        end do

    end if

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
!The interpolation used here is quadrilinear interpolation which is a straightforward
!generalization of bilinear interpolation.
!
!The function also allows extrapolation and does not perform any checks
!as to whether we are between grid boundaries (to increase efficiency). We should bever
!use this subroutine for extrapolation, though, because that leads to all sorts
!of problems.
!It uses a guess of position which is suplied in j_guess. The point is that
!in maximization, we usually evaluate points which are close to each other,
!so we can have a good idea of where the neighbouring points in the grid are.
!This can save a lot of time (my estimate is that with 100 grid points per variable,
!the time to find the point should be cut by a factor of about 81), so it is
!worth the effort (interpolation is very deep in the program, we need to be as
!efficient as possible).
subroutine interp_V(V,x,j,s_ind,y,a1_gr,a2_gr,rho_gr,aex_gr,N_a,N_rho,N_aex)
    real(dp), dimension(:,:,:,:,:), intent(in) :: V
    real(dp), dimension(4), intent(in) :: x
    integer, dimension(4), intent(inout) :: j !on entry, the guess of position of
    !the point x on the grid. On exit, the actual position. If the first element
    !if zero (j(1) == 0), then we do not use the initial guess, and use subroutine
    !locate_gr rather than hunt (this subroutine is more efficient if we do not have a good
    !initial guess).

    integer, intent(in) :: s_ind
    real(dp), intent(out) :: y
    real(dp), dimension(:) ,intent(in) :: a1_gr,a2_gr,rho_gr,aex_gr

    integer, intent(in) :: N_a,N_rho,N_aex
    !Same number of gridpoints for agent 1 and agent 2 is assumed - this can be
    !generalized easily, but its difficult to see the benefit of doing it.

    real(dp) :: t_a1,t_a2,t_rho,t_aex ! weights in linear combinations t*x(j) + (1-t)*x(j+1)

    !This follows quite closely algorithms suggested in the Numerical recipes in Fortran book.
    !1) Find the points to be used in interpolation (for all variables, point with index j
    !is the point s.t. x lies between x(j) and x(j+1).

    if(j(1) == 0) then
        call locate_gr(a1_gr,N_a,x(1),j(1)) !j_a1 will be the index of the grid point in a1_gr
        call locate_gr(a2_gr,N_a,x(2),j(2))
        call locate_gr(rho_gr,N_rho,x(3),j(3))
        call locate_gr(aex_gr,N_aex,x(4),j(4))
    else
        call hunt(a1_gr,N_a,x(1),j(1))
        call hunt(a2_gr,N_a,x(2),j(2))
        call hunt(rho_gr,N_rho,x(3),j(3))
        call hunt(aex_gr,N_aex,x(4),j(4))
    end if

    !If we are at the border of the grid in any dimension, the returned index might be 0 (it's
    !down to numerical errors). Then we would have a problem (trying to evaluate V(0,1,1,1,1) for example
    !The program would not crash, but it would be some non-sensical value - usually value at a different
    !point in state space. Therefore, floor all the j's at value 1.
    j(1) = max(1,j(1))
    j(2) = max(1,j(2))
    j(3) = max(1,j(3))
    j(4) = max(1,j(4))

    !coefficients in linear combinations (value t = 0 is associated with the higher point (j+1))
    t_a1 = (a1_gr(j(1)+1)-x(1))/(a1_gr(j(1)+1)-a1_gr(j(1)))
    t_a2 = (a2_gr(j(2)+1)-x(2))/(a2_gr(j(2)+1)-a2_gr(j(2)))
    t_rho = (rho_gr(j(3)+1)-x(3))/(rho_gr(j(3)+1)-rho_gr(j(3)))
    t_aex = (aex_gr(j(4)+1)-x(4))/(aex_gr(j(4)+1)-aex_gr(j(4)))

    !The following looks horrible but it should be ok (generated using a spreadsheet script)
    !(this subroutine is used potentially millions of times per one iteration in VFI algorithm
    !so maximum efficiency is important - no unnecesary saving of intermediate results or calling
    !subroutines). Also, elements of V are accessed in the order in which they are stored
    !in memory - it could be the case that the compiler does optimization of this sort automatically
    !but it doesn't hurt to write it explicitly.
    y = V(j(1),j(2),j(3),j(4),s_ind)*t_a1*t_a2*t_rho*t_aex+&
        V(j(1)+1,j(2),j(3),j(4),s_ind)*(1-t_a1)*t_a2*t_rho*t_aex+&
        V(j(1),j(2)+1,j(3),j(4),s_ind)*t_a1*(1-t_a2)*t_rho*t_aex+&
        V(j(1)+1,j(2)+1,j(3),j(4),s_ind)*(1-t_a1)*(1-t_a2)*t_rho*t_aex+&
        V(j(1),j(2),j(3)+1,j(4),s_ind)*t_a1*t_a2*(1-t_rho)*t_aex+&
        V(j(1)+1,j(2),j(3)+1,j(4),s_ind)*(1-t_a1)*t_a2*(1-t_rho)*t_aex+&
        V(j(1),j(2)+1,j(3)+1,j(4),s_ind)*t_a1*(1-t_a2)*(1-t_rho)*t_aex+&
        V(j(1)+1,j(2)+1,j(3)+1,j(4),s_ind)*(1-t_a1)*(1-t_a2)*(1-t_rho)*t_aex+&
        V(j(1),j(2),j(3),j(4)+1,s_ind)*t_a1*t_a2*t_rho*(1-t_aex)+&
        V(j(1)+1,j(2),j(3),j(4)+1,s_ind)*(1-t_a1)*t_a2*t_rho*(1-t_aex)+&
        V(j(1),j(2)+1,j(3),j(4)+1,s_ind)*t_a1*(1-t_a2)*t_rho*(1-t_aex)+&
        V(j(1)+1,j(2)+1,j(3),j(4)+1,s_ind)*(1-t_a1)*(1-t_a2)*t_rho*(1-t_aex)+&
        V(j(1),j(2),j(3)+1,j(4)+1,s_ind)*t_a1*t_a2*(1-t_rho)*(1-t_aex)+&
        V(j(1)+1,j(2),j(3)+1,j(4)+1,s_ind)*(1-t_a1)*t_a2*(1-t_rho)*(1-t_aex)+&
        V(j(1),j(2)+1,j(3)+1,j(4)+1,s_ind)*t_a1*(1-t_a2)*(1-t_rho)*(1-t_aex)+&
        V(j(1)+1,j(2)+1,j(3)+1,j(4)+1,s_ind)*(1-t_a1)*(1-t_a2)*(1-t_rho)*(1-t_aex)
end subroutine interp_V

!Subroutine interp_V_shep is just like interp_V, but a local Shepard's interpolation
!is used. There is one extra parameter shep_n, which tells the program how many
!closest points in every direction on the grid are to be used. For example, if
!shep_n = 2, then only 2^4 points are used (but it makes pretty much zero sense to
!use this if shep_n = 2, because it is slower than quadrilinear inrepolation, and
!there is no good reason to expect this to perform better). Reasonable values
!would be between 3 and 5 (larger values lead to slow interpolation).
!
!The weights used are the same as in NaG subroutine for Shepard's N-dimensional
!interpolation. I coded this manually because that subroutine does not assume
!an ordered grid and is less efficient in this particular case.
subroutine interp_V_shep(V,x,j,s_ind,y,a1_gr,a2_gr,rho_gr,aex_gr,N_a,N_rho,N_aex,shep_n,shep_norm,shep_norm_val,shep_floor)
    real(dp), dimension(:,:,:,:,:), intent(in) :: V
    real(dp), dimension(4), intent(in) :: x
    integer, dimension(4), intent(inout) :: j !on entry, the guess of position of
    !the point x on the grid. On exit, the actual position. If the first element
    !if zero (j(1) == 0), then we do not use the initial guess, and use subroutine
    !locate_gr rather than hunt (this subroutine is more efficient if we do not have a good
    !initial guess).

    integer, intent(in) :: s_ind
    real(dp), intent(out) :: y
    real(dp), dimension(:) ,intent(in) :: a1_gr,a2_gr,rho_gr,aex_gr

    integer, intent(in) :: N_a,N_rho,N_aex
    !Same number of gridpoints for agent 1 and agent 2 is assumed - this can be
    !generalized easily, but its difficult to see the benefit of doing it.
    integer, intent(in) :: shep_n !number of grid points to be used in every dimension
    logical, intent(in) :: shep_norm !if this is tru, then the distance in each dimension
    !is taken as share of overall distance covered by the grid (this way, points are
    !more penalized (in the sense of lower weight) if they are distant in rho dimension
    !(whose grid is narrow and which varies only little), than if they are equally
    !distant in a different dimension.
    real(dp), dimension(4), intent(in) :: shep_norm_val
    real(dp), intent(in) :: shep_floor

    integer :: j_min(4) !For each of the 4 variables, this is the index of the
    !first element to be used in the interpolation

    integer, dimension(:,:), pointer :: shep_inc_pntr
    integer :: inc_ind !index of increment
    integer, dimension(4) :: ind_all !all indices (used repeatedly so it's worth saving it)

    real(dp) :: sow !sum of weights
    real(dp) :: weight

    !The pre-computed matrix of increments (pointer associated in union.f90)
    common /shep_pointer/ shep_inc_pntr

    !This follows quite closely algorithms suggested in the Numerical recipes in Fortran book.
    !(1) Find the points to be used in interpolation (for all variables, point with index j
    !is the point s.t. x lies between x(j) and x(j+1).

    if(j(1) == 0) then
        call locate_gr(a1_gr,N_a,x(1),j(1)) !j_a1 will be the index of the grid point in a1_gr
        call locate_gr(a2_gr,N_a,x(2),j(2))
        call locate_gr(rho_gr,N_rho,x(3),j(3))
        call locate_gr(aex_gr,N_aex,x(4),j(4))
    else
        call hunt(a1_gr,N_a,x(1),j(1))
        call hunt(a2_gr,N_a,x(2),j(2))
        call hunt(rho_gr,N_rho,x(3),j(3))
        call hunt(aex_gr,N_aex,x(4),j(4))
    end if

    !If we are at the border of the grid in any dimension, the returned index might be 0 (it's
    !down to numerical errors). Then we would have a problem (trying to evaluate V(0,1,1,1,1) for example
    !The program would not crash, but it would be some non-sensical value - usually value at a different
    !point in state space. Therefore, floor all the j's at value 1.
    j(1) = max(1,j(1))
    j(2) = max(1,j(2))
    j(3) = max(1,j(3))
    j(4) = max(1,j(4))

    !(2) Now (for each dimension) find the index of the first point in the grid to be used
    !in interpolation (this is truncated by the bounds of the grid indices).
    j_min = min(max(j-(shep_n-1)/2,1),[N_a,N_a,N_rho,N_aex]+1-shep_n)

    !(3) For each used point, compute the weight, its contribution to sum of weights,
    !and it's contribution to the value

    !I do this with one cycle only for extra efficiency
    !For the first variable (a1) the range is j_min(1) + shep_n - 1. Similarly for the other variables

    !Initialize sum of weights and the interpolated value
    y = 0.0_dp
    sow = 0.0_dp

    !IF shep_norm == .true., we normalize the distance in every dimension by the
    !total range of grid in that dimension

    if(shep_norm) then
        do inc_ind = 1,shep_n**4
            !The indices of the points in terms of the value function are
            !j_min + shep_inc_pntr(inc_ind,:). We use this repeatedly so it is worth saving
            ind_all = j_min + shep_inc_pntr(inc_ind,:)

            !A sensible way to do this would be to scale each of the deviations
            !by the total range of bounds for that variable. The current version does not do this.

            weight = 1.0_dp/(sum((([x(1),x(2),x(3),x(4)] - &
            [a1_gr(ind_all(1)),a2_gr(ind_all(2)),rho_gr(ind_all(3)),aex_gr(ind_all(4))])/shep_norm_val)&
                ** 2.0_dp) + shep_floor)
                 !(we added the shep_floor so that we don't divide by zero when the point happens to fall
                 !exactly on the grid)

            !Add the contribution of the point to the interpolated value
            y = y + weight * V(ind_all(1),ind_all(2),ind_all(3),ind_all(4),s_ind)
            sow = sow + weight
        end do
    else

        do inc_ind = 1,shep_n**4
            !The indices of the points in terms of the value function are
            !j_min + shep_inc_pntr(inc_ind,:). We use this repeatedly so it is worth saving
            ind_all = j_min + shep_inc_pntr(inc_ind,:)

            weight = 1.0_dp/(sum(([x(1),x(2),x(3),x(4)] - [a1_gr(ind_all(1)),a2_gr(ind_all(2)),rho_gr(ind_all(3)),aex_gr(ind_all(4))])&
                ** 2.0_dp) + shep_floor)

            !Add the contribution of the point to the interpolated value
            y = y + weight * V(ind_all(1),ind_all(2),ind_all(3),ind_all(4),s_ind)
            sow = sow + weight
        end do

    end if
    y = y / sow !normalize (so that the weights sum to 1)

end subroutine interp_V_shep

subroutine interp_V_spline(V,x,s_ind,y,a1_gr,a2_gr,rho_gr,aex_gr,N_a,N_rho,N_aex,&
bspline_k_a1,bspline_k_a2,bspline_k_rho,bspline_k_aex)
    real(dp), dimension(:,:,:,:,:), intent(in) :: V
    real(dp), dimension(4), intent(in) :: x
    integer, intent(in) :: s_ind
    real(dp), intent(out) :: y
    real(dp), dimension(:) ,intent(in) :: a1_gr,a2_gr,rho_gr,aex_gr
    integer, intent(in) :: N_a,N_rho,N_aex,bspline_k_a1,bspline_k_a2,bspline_k_rho,bspline_k_aex

    real(dp), dimension(:,:), pointer :: t_a1_pntr,t_a2_pntr,t_rho_pntr,t_aex_pntr
    integer, dimension(M_par) :: inbvx,inbvy,inbvz,inbvq,iloy,iloz,iloq
    !initialization parameter which must be set to 1 the first time this routine is called
    !(every time that I recompute the bspline coefficients).

    integer :: bspline_iflag !if this is not zero then some error happened.

    !Some inputs are passed to this subroutine directly from the main program using a common block.
    common /bspline_pointer1/ inbvx, inbvy, inbvz, inbvq, iloy, iloz, iloq,&
    t_a1_pntr,t_a2_pntr,t_rho_pntr,t_aex_pntr


call db4val(x(1),x(2),x(3),x(4),&
                                0,0,0,0,& !the zeroes here mean that we want the interpolant itself (0th derivative)
                                t_a1_pntr(1:(N_a + bspline_k_a1),s_ind),t_a2_pntr(1:(N_a + bspline_k_a2),s_ind)&
                                ,t_rho_pntr(1:(N_a + bspline_k_rho),s_ind),t_aex_pntr(1:(N_aex + bspline_k_aex),s_ind),&
                                N_a,N_a,N_rho,N_aex,&
                                bspline_k_a1,bspline_k_a2,bspline_k_rho,bspline_k_aex,&
                                V(:,:,:,:,s_ind),y,bspline_iflag,&
                                inbvx(s_ind),inbvy(s_ind),inbvz(s_ind),inbvq(s_ind),&
                                iloy(s_ind),iloz(s_ind),iloq(s_ind),.true.)
!(I leave the last argument as .true. (allowing extrapolation) - this should never be used in practice
!but it might prevent some crashes due to rounding errors).

end subroutine interp_V_spline


subroutine interp_x_shep(x_pol_all_unf,x,j,s_ind,y,a1_gr,a2_gr,rho_gr,aex_gr,N_a,N_rho,N_aex,shep_n,shep_norm,shep_norm_val&
,shep_floor)
    real(dp), dimension(:,:,:,:,:,:), intent(in) :: x_pol_all_unf !policy function in 'unfolded' form. First
    !index corresponds to value of choice variables (there are ncvar_par of them), the following
    !5 indices tell us the index of the point at which the value of policy function contained in the
    !first element is attained.
    real(dp), dimension(4), intent(in) :: x !(a1,a2,rho,aex) - current choice

    integer, dimension(4), intent(inout) :: j !on entry, the guess of position of
    !the point x on the grid. On exit, the actual position. If the first element
    !if zero (j(1) == 0), then we do not use the initial guess, and use subroutine
    !locate_gr rather than hunt (this subroutine is more efficient if we do not have a good
    !initial guess).

    integer, intent(in) :: s_ind
    real(dp), intent(out), dimension(ncvar_par) :: y !interpolated choice vector
    real(dp), dimension(:), intent(in) :: a1_gr,a2_gr,rho_gr,aex_gr
    real(dp), intent(in) :: shep_floor

    integer, intent(in) :: N_a,N_rho,N_aex
    !Same number of gridpoints for agent 1 and agent 2 is assumed - this can be
    !generalized easily, but its difficult to see the benefit of doing it.

    integer, intent(in) :: shep_n !number of grid points to be used in every dimension
    logical, intent(in) :: shep_norm !if this is tru, then the distance in each dimension
    !is taken as share of overall distance covered by the grid (this way, points are
    !more penalized (in the sense of lower weight) if they are distant in rho dimension
    !(whose grid is narrow and which varies only little), than if they are equally
    !distant in a different dimension.
    real(dp), dimension(4), intent(in) :: shep_norm_val

    integer :: j_min(4) !For each of the 4 variables, this is the index of the
    !first element to be used in the interpolation


    integer, dimension(:,:), pointer :: shep_inc_pntr
    integer :: inc_ind !index of increment
    integer, dimension(4) :: ind_all !all indices (used repeatedly so it's worth saving it)

    real(dp) :: sow !sum of weights
    real(dp) :: weight

    !The pre-computed matrix of increments (pointer associated in union.f90)
    common /shep_pointer/ shep_inc_pntr

    integer :: m_ind,i_ind !indices for cycling over elements of choice matrix.

    !This follows quite closely algorithm suggested in the Numerical recipes in Fortran book.
    !1) Find the points to be used in interpolation (for all variables, point with index j
    !is the point s.t. x lies between x(j) and x(j+1).
    if(j(1) == 0) then
        call locate_gr(a1_gr,N_a,x(1),j(1)) !j_a1 will be the index of the grid point in a1_gr
        call locate_gr(a2_gr,N_a,x(2),j(2))
        call locate_gr(rho_gr,N_rho,x(3),j(3))
        call locate_gr(aex_gr,N_aex,x(4),j(4))
    else
        call hunt(a1_gr,N_a,x(1),j(1))
        call hunt(a2_gr,N_a,x(2),j(2))
        call hunt(rho_gr,N_rho,x(3),j(3))
        call hunt(aex_gr,N_aex,x(4),j(4))
    end if

    !If we are at the border of the grid in any dimension, the returned index might be 0 (it's
    !down to numerical errors). Then we would have a problem (trying to evaluate V(0,1,1,1,1) for example
    !The program would not crash, but it would be some non-sensical value - usually value at a different
    !point in state space. Therefore, floor all the j's at value 1.
    j(1) = max(1,j(1))
    j(2) = max(1,j(2))
    j(3) = max(1,j(3))
    j(4) = max(1,j(4))

    !(2) Now (for each dimension) find the index of the first point in the grid to be used
    !in interpolation (this is truncated by the bounds of the grid indices).
    j_min = min(max(j-(shep_n-1)/2,1),[N_a,N_a,N_rho,N_aex]+1-shep_n)

    !(3) For each used point, compute the weight, its contribution to sum of weights,
    !and it's contribution to the value

    !I do this with one cycle only for extra efficiency
    !For the first variable (a1) the range is j_min(1) + shep_n - 1. Similarly for the other variables

    !Initialize sum of weights and the interpolated value
    y = 0.0_dp
    sow = 0.0_dp

    if(shep_norm) then
        do inc_ind = 1,shep_n**4
            !The indices of the points in terms of the value function are
            !j_min + shep_inc_pntr(inc_ind,:). We use this repeatedly so it is worth saving
            ind_all = j_min + shep_inc_pntr(inc_ind,:)

            !A sensible way to do this would be to scale each of the deviations
            !by the total range of bounds for that variable. The current version does not do this.

            weight = 1.0_dp/(sum((([x(1),x(2),x(3),x(4)] - &
            [a1_gr(ind_all(1)),a2_gr(ind_all(2)),rho_gr(ind_all(3)),aex_gr(ind_all(4))])/shep_norm_val)&
                ** 2.0_dp) + shep_floor)
                 !(we added the shep_floor so that we don't divide by zero when the point happens to fall
                 !exactly on the grid).


            !Add the contribution of the point to the interpolated policy
            y = y + weight*x_pol_all_unf(:,ind_all(1),ind_all(2),ind_all(3),ind_all(4),s_ind)

            sow = sow + weight
        end do
    else

        do inc_ind = 1,shep_n**4
            !The indices of the points in terms of the value function are
            !j_min + shep_inc_pntr(inc_ind,:). We use this repeatedly so it is worth saving
            ind_all = j_min + shep_inc_pntr(inc_ind,:)

            weight = 1.0_dp/(sum(([x(1),x(2),x(3),x(4)] - [a1_gr(ind_all(1)),a2_gr(ind_all(2)),rho_gr(ind_all(3)),aex_gr(ind_all(4))])&
                ** 2.0_dp) + shep_floor)
                 !(we added the 10E-12 so that we don't divide by zero when the point happens to fall
                 !exactly on the grid)

            !Add the contribution of the point to the interpolated value
            y = y + weight*x_pol_all_unf(:,ind_all(1),ind_all(2),ind_all(3),ind_all(4),s_ind)
            sow = sow + weight
        end do

    end if
    y = y / sow !normalize (so that the weights sum to 1)

end subroutine interp_x_shep

!A debug version of interp_x_shep which writes a lot of additional information.
subroutine interp_x_shep_debug(pars,x_pol_all_unf,x,j,s_ind,y,a1_gr,a2_gr,rho_gr,aex_gr,N_a,N_rho,N_aex,shep_n,shep_norm,shep_norm_val&
,shep_floor)
    real(dp), dimension(:,:,:,:,:,:), intent(in) :: x_pol_all_unf !policy function in 'unfolded' form. First
    !index corresponds to value of choice variables (there are ncvar_par of them), the following
    !5 indices tell us the index of the point at which the value of policy function contained in the
    !first element is attained.
    real(dp), dimension(4), intent(in) :: x !(a1,a2,rho,aex) - current choice

    type(par), intent(in) :: pars

    integer, dimension(4), intent(inout) :: j !on entry, the guess of position of
    !the point x on the grid. On exit, the actual position. If the first element
    !if zero (j(1) == 0), then we do not use the initial guess, and use subroutine
    !locate_gr rather than hunt (this subroutine is more efficient if we do not have a good
    !initial guess).

    integer, intent(in) :: s_ind
    real(dp), intent(out), dimension(ncvar_par) :: y !interpolated choice vector
    real(dp), dimension(:), intent(in) :: a1_gr,a2_gr,rho_gr,aex_gr

    integer, intent(in) :: N_a,N_rho,N_aex
    !Same number of gridpoints for agent 1 and agent 2 is assumed - this can be
    !generalized easily, but its difficult to see the benefit of doing it.

    integer, intent(in) :: shep_n !number of grid points to be used in every dimension
    logical, intent(in) :: shep_norm !if this is tru, then the distance in each dimension
    !is taken as share of overall distance covered by the grid (this way, points are
    !more penalized (in the sense of lower weight) if they are distant in rho dimension
    !(whose grid is narrow and which varies only little), than if they are equally
    !distant in a different dimension.
    real(dp), dimension(4), intent(in) :: shep_norm_val
    real(dp), intent(in) :: shep_floor

    integer :: j_min(4) !For each of the 4 variables, this is the index of the
    !first element to be used in the interpolation

    !tmp debug variable
    real(dp), dimension(ncvaR_par) :: xc


    integer, dimension(:,:), pointer :: shep_inc_pntr
    integer :: inc_ind !index of increment
    integer, dimension(4) :: ind_all !all indices (used repeatedly so it's worth saving it)

    real(dp) :: sow !sum of weights
    real(dp) :: weight


    real(dp), dimension(I_par) :: a_st
    real(dp), dimension(I_par-1) :: rho_st
    real(dp) :: aex_st
    integer :: g_st_ind

    !Current_period choices (usually state-contingent plans. Only one element of these is saved into sim_series,
    !according to shock realization).
    real(dp), dimension(M_par,I_par) :: u_c,u_l
    real(dp), dimension(M_par,1) :: R
    real(dp), dimension(M_par,I_par) :: c,l !consumption and labour supply
    real(dp), dimension(M_par,I_par) :: a_pr !next-period adjusted MU-asset holdings
    real(dp), dimension(M_par,1) :: rho_pr !next-period ratio of marginal utilities
    real(dp), dimension(M_par,1) :: aex_pr !next-period MU-adjusted external debt of country 1
    integer :: getavfail

    real(dp), dimension(M_par,I_par) :: tau_plan


    !The pre-computed matrix of increments (pointer associated in union.f90)
    common /shep_pointer/ shep_inc_pntr

    integer :: m_ind,i_ind !indices for cycling over elements of choice matrix.

    !This follows quite closely algorithm suggested in the Numerical recipes in Fortran book.
    !1) Find the points to be used in interpolation (for all variables, point with index j
    !is the point s.t. x lies between x(j) and x(j+1).
    if(j(1) == 0) then
        call locate_gr(a1_gr,N_a,x(1),j(1)) !j_a1 will be the index of the grid point in a1_gr
        call locate_gr(a2_gr,N_a,x(2),j(2))
        call locate_gr(rho_gr,N_rho,x(3),j(3))
        call locate_gr(aex_gr,N_aex,x(4),j(4))
    else
        call hunt(a1_gr,N_a,x(1),j(1))
        call hunt(a2_gr,N_a,x(2),j(2))
        call hunt(rho_gr,N_rho,x(3),j(3))
        call hunt(aex_gr,N_aex,x(4),j(4))
    end if

    !If we are at the border of the grid in any dimension, the returned index might be 0 (it's
    !down to numerical errors). Then we would have a problem (trying to evaluate V(0,1,1,1,1) for example
    !The program would not crash, but it would be some non-sensical value - usually value at a different
    !point in state space. Therefore, floor all the j's at value 1.
    j(1) = max(1,j(1))
    j(2) = max(1,j(2))
    j(3) = max(1,j(3))
    j(4) = max(1,j(4))

    !(2) Now (for each dimension) find the index of the first point in the grid to be used
    !in interpolation (this is truncated by the bounds of the grid indices).
    j_min = min(max(j-(shep_n-1)/2,1),[N_a,N_a,N_rho,N_aex]+1-shep_n)

    !(3) For each used point, compute the weight, its contribution to sum of weights,
    !and it's contribution to the value

    !I do this with one cycle only for extra efficiency
    !For the first variable (a1) the range is j_min(1) + shep_n - 1. Similarly for the other variables

    !Initialize sum of weights and the interpolated value
    y = 0.0_dp
    sow = 0.0_dp

    write(*,*) '=============debug in interp_x_shep_debug================='

    if(shep_norm) then
        do inc_ind = 1,shep_n**4
            !The indices of the points in terms of the value function are
            !j_min + shep_inc_pntr(inc_ind,:). We use this repeatedly so it is worth saving
            ind_all = j_min + shep_inc_pntr(inc_ind,:)

            !A sensible way to do this would be to scale each of the deviations
            !by the total range of bounds for that variable. The current version does not do this.

            weight = 1.0_dp/(sum((([x(1),x(2),x(3),x(4)] - &
            [a1_gr(ind_all(1)),a2_gr(ind_all(2)),rho_gr(ind_all(3)),aex_gr(ind_all(4))])/shep_norm_val)&
                ** 2.0_dp) + shep_floor)
                 !(we added the 10E-12 so that we don't divide by zero when the point happens to fall
                 !exactly on the grid)

            !Add the contribution of the point to the interpolated policy
            y = y + weight*x_pol_all_unf(:,ind_all(1),ind_all(2),ind_all(3),ind_all(4),s_ind)

            sow = sow + weight
        end do
    else

        do inc_ind = 1,shep_n**4
            !The indices of the points in terms of the value function are
            !j_min + shep_inc_pntr(inc_ind,:). We use this repeatedly so it is worth saving
            ind_all = j_min + shep_inc_pntr(inc_ind,:)

            weight = 1.0_dp/(sum(([x(1),x(2),x(3),x(4)] - [a1_gr(ind_all(1)),a2_gr(ind_all(2)),rho_gr(ind_all(3)),aex_gr(ind_all(4))])&
                ** 2.0_dp) + shep_floor)
                 !(we added the 10E-12 so that we don't divide by zero when the point happens to fall
                 !exactly on the grid)

            !Add the contribution of the point to the interpolated value
            y = y + weight*x_pol_all_unf(:,ind_all(1),ind_all(2),ind_all(3),ind_all(4),s_ind)
            sow = sow + weight

            !Now we write the value of policy function at all the points in question.
            !this already assumes that we are running the program on 1 CPU only.

            write(*,*) 'point number ',inc_ind
            write(*,*) 'value of states (a1,a2,rho,aex):'
            write(*,*) [a1_gr(ind_all(1)),a2_gr(ind_all(2)),rho_gr(ind_all(3)),aex_gr(ind_all(4))]
            write(*,*) 's_ind = ',s_ind
            write(*,*) 'Policy at this point:'
            write(*,*) x_pol_all_unf(:,ind_all(1),ind_all(2),ind_all(3),ind_all(4),s_ind)

            xc = x_pol_all_unf(:,ind_all(1),ind_all(2),ind_all(3),ind_all(4),s_ind)

            !Now get all the remaining variables.

            a_st(1) = a1_gr(ind_all(1))
            a_st(2) = a2_gr(ind_all(2))
            rho_st = rho_gr(ind_all(3))
            aex_st = aex_gr(ind_all(4))

            !call a debug version of get_av (which write a lot of output).
            call get_av(pars,xc,a_st,rho_st,aex_st,s_ind,getavfail,c,l,a_pr,rho_pr,aex_pr,.false.)

            !To get more debug info just use get_av_debug

            call util_c(c,u_c,pars%A_par,pars%sigma_par)
            call util_l(l,u_l,pars%B_par,pars%Gamma_par)

            !Now get the tax rate in both countries for both shock realizations.

            !later add here the value of tax rate in both countries - I will need to call get_av.
        !here assume that productivity is 1 in both countries
        tau_plan = u_l/u_c + 1.0_dp

        write(*,*) 'tau (state 1) = ',tau_plan(1,:)
        write(*,*) 'tau (state 2) = ',tau_plan(2,:)
        write(*,*) 'tau (state 1 - state 2) = ', tau_plan(1,:) - tau_plan(2,:)

!        write(*,*) 'c = ',c
!        write(*,*) 'l = ',l
        write(*,*) 'rho_pr = ',rho_pr

        write(*,*) 'weight = ',weight

!        write(*,*) 'stopping after the first point'
!        error stop



            write(*,*) '________________________________________________'


        end do

    end if
    y = y / sow !normalize (so that the weights sum to 1)

end subroutine interp_x_shep_debug

!subroutine interp_x serves for interpolating policy function - it is essentially pretty much the same
!as interp_V, the difference being that we get an interpolated vector rather than a number.
!The interpolation is done element by element.
!
!This subroutine does not use an initial guess as interp_V, and is thus less efficient.
!We can improve this at some point but it's a low priority, because interpolation of policy
!function is used very rarely (mainly in between stages of CTS algorithm), and it is not
!clear that we could find a good initial guess in this case anyway.
subroutine interp_x(x_pol_all_unf,x,j,s_ind,y,a1_gr,a2_gr,rho_gr,aex_gr,N_a,N_rho,N_aex)
    real(dp), dimension(:,:,:,:,:,:), intent(in) :: x_pol_all_unf !policy function in 'unfolded' form. First
    !index corresponds to value of choice variables (there are ncvar_par of them), the following
    !5 indices tell us the index of the point at which the value of policy function contained in the
    !first element is attained.
    real(dp), dimension(4), intent(in) :: x !(a1,a2,rho,aex) - current choice
    integer, intent(in) :: s_ind
    real(dp), intent(out), dimension(ncvar_par) :: y !interpolated choice vector
    real(dp), dimension(:), intent(in) :: a1_gr,a2_gr,rho_gr,aex_gr

    integer, dimension(4), intent(inout) :: j !on entry, the guess of position of
    !the point x on the grid. On exit, the actual position. If the first element
    !if zero (j(1) == 0), then we do not use the initial guess, and use subroutine
    !locate_gr rather than hunt (this subroutine is more efficient if we do not have a good
    !initial guess).

    integer, intent(in) :: N_a,N_rho,N_aex
    !Same number of gridpoints for agent 1 and agent 2 is assumed - this can be
    !generalized easily, but its difficult to see the benefit of doing it.
    integer :: j_a1,j_a2,j_rho,j_aex !used for storing which points are to be used in interpolation.
    real(dp) :: t_a1,t_a2,t_rho,t_aex ! weights in linear combinations t*x(j) + (1-t)*x(j+1)

    integer :: m_ind,i_ind !indices for cycling over elements of choice matrix.



    y = 0.0_dp !initilalize

    !1) Find the points to be used in interpolation (for all variables, point with index j
    !is the point s.t. x lies between x(j) and x(j+1).
    if(j(1) == 0) then
        call locate_gr(a1_gr,N_a,x(1),j(1)) !j_a1 will be the index of the grid point in a1_gr
        call locate_gr(a2_gr,N_a,x(2),j(2))
        call locate_gr(rho_gr,N_rho,x(3),j(3))
        call locate_gr(aex_gr,N_aex,x(4),j(4))
    else
        call hunt(a1_gr,N_a,x(1),j(1))
        call hunt(a2_gr,N_a,x(2),j(2))
        call hunt(rho_gr,N_rho,x(3),j(3))
        call hunt(aex_gr,N_aex,x(4),j(4))
    end if


    j(1) = max(1,j(1))
    j(2) = max(1,j(2))
    j(3) = max(1,j(3))
    j(4) = max(1,j(4))


    !coefficients in linear combinations (value 0 is associated with the lower point)
    t_a1 = (a1_gr(j(1)+1)-x(1))/(a1_gr(j(1)+1)-a1_gr(j(1)))
    t_a2 = (a2_gr(j(2)+1)-x(2))/(a2_gr(j(2)+1)-a2_gr(j(2)))
    t_rho = (rho_gr(j(3)+1)-x(3))/(rho_gr(j(3)+1)-rho_gr(j(3)))
    t_aex = (aex_gr(j(4)+1)-x(4))/(aex_gr(j(4)+1)-aex_gr(j(4)))

    !The interpolation for each of the elements of the consumption vector is done independently.
        y = x_pol_all_unf(:,j(1),j(2),j(3),j(4),s_ind)*t_a1*t_a2*t_rho*t_aex+&
            x_pol_all_unf(:,j(1)+1,j(2),j(3),j(4),s_ind)*(1-t_a1)*t_a2*t_rho*t_aex+&
            x_pol_all_unf(:,j(1),j(2)+1,j(3),j(4),s_ind)*t_a1*(1-t_a2)*t_rho*t_aex+&
            x_pol_all_unf(:,j(1)+1,j(2)+1,j(3),j(4),s_ind)*(1-t_a1)*(1-t_a2)*t_rho*t_aex+&
            x_pol_all_unf(:,j(1),j(2),j(3)+1,j(4),s_ind)*t_a1*t_a2*(1-t_rho)*t_aex+&
            x_pol_all_unf(:,j(1)+1,j(2),j(3)+1,j(4),s_ind)*(1-t_a1)*t_a2*(1-t_rho)*t_aex+&
            x_pol_all_unf(:,j(1),j(2)+1,j(3)+1,j(4),s_ind)*t_a1*(1-t_a2)*(1-t_rho)*t_aex+&
            x_pol_all_unf(:,j(1)+1,j(2)+1,j(3)+1,j(4),s_ind)*(1-t_a1)*(1-t_a2)*(1-t_rho)*t_aex+&
            x_pol_all_unf(:,j(1),j(2),j(3),j(4)+1,s_ind)*t_a1*t_a2*t_rho*(1-t_aex)+&
            x_pol_all_unf(:,j(1)+1,j(2),j(3),j(4)+1,s_ind)*(1-t_a1)*t_a2*t_rho*(1-t_aex)+&
            x_pol_all_unf(:,j(1),j(2)+1,j(3),j(4)+1,s_ind)*t_a1*(1-t_a2)*t_rho*(1-t_aex)+&
            x_pol_all_unf(:,j(1)+1,j(2)+1,j(3),j(4)+1,s_ind)*(1-t_a1)*(1-t_a2)*t_rho*(1-t_aex)+&
            x_pol_all_unf(:,j(1),j(2),j(3)+1,j(4)+1,s_ind)*t_a1*t_a2*(1-t_rho)*(1-t_aex)+&
            x_pol_all_unf(:,j(1)+1,j(2),j(3)+1,j(4)+1,s_ind)*(1-t_a1)*t_a2*(1-t_rho)*(1-t_aex)+&
            x_pol_all_unf(:,j(1),j(2)+1,j(3)+1,j(4)+1,s_ind)*t_a1*(1-t_a2)*(1-t_rho)*(1-t_aex)+&
            x_pol_all_unf(:,j(1)+1,j(2)+1,j(3)+1,j(4)+1,s_ind)*(1-t_a1)*(1-t_a2)*(1-t_rho)*(1-t_aex)
end subroutine interp_x




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

!Subroutine hunt is a subroutine for locating position of a point in a grid of ascending
!order, using an initial guess of position (in this regard, it is more efficient than
!subroutine locate_gr). This is copied from Numerical Recipes in Fortran 77 book, with minor
!modifications. For example, we do not check that the order of the grid is ascending
!(we assume that the grids have been constructed properly).
subroutine hunt(xx,n,x,jlo)
integer :: jlo,n
real(dp) ::  x,xx(n)
!Given an array xx(1:n), and given a value x , returns a value jlo such that x is between
!xx(jlo) and xx(jlo+1) . xx(1:n) must be monotonic, either increasing or decreasing.
!jlo=0 or jlo=n is returned to indicate that x is out of range. jlo on input is taken as
!the initial guess for jlo on output.
integer :: inc,jhi,jm

!This checking of whether the initial guess is in the bounds should always be unnecessary
!(but it's safer to keep it here).
if(jlo <= 0 .or. jlo > n) then
    !Input guess not useful. Go immediately to bisection.
    jlo=0
    jhi=n+1
    goto 3
endif

inc=1 !Set the hunting increment (initially = 1).

if(x >= xx(jlo)) then
    !Hunt up:
1   jhi=jlo+inc
    if(jhi > n)then
        !Done hunting, since off end of table.
        jhi=n+1
    else if(x >= xx(jhi))then
        !Not done hunting, so double the increment and try again.
        jlo=jhi
        inc=inc+inc
        goto 1
    endif
!Done hunting, value bracketed.
else
    !Hunt down:
    jhi=jlo
2   jlo=jhi-inc
    if(jlo < 1)then
        !Done hunting, since off end of table.
        jlo=0
    else if(x < xx(jlo))then
        !Not done hunting, so double the increment and try again.
        jhi=jlo
        inc=inc+inc
        goto 2
    endif
!Done hunting, value bracketed.
endif
!Hunt is done, so begin the final bisection phase:
3 if(jhi - jlo == 1)then
    if(x >= xx(n)) jlo=n-1
    if(x <= xx(1)) jlo=1
    return
endif
jm=(jhi+jlo)/2
if(x >= xx(jm))then
    jlo=jm
else
    jhi=jm
endif
goto 3
end subroutine hunt


!subroutine LFFC (look for feasible choice) attempts to find a choice of consumption
!such that no constraint is violated, including the constraints imposed by grid range
!of a and rho. It does so for a given value of states (in the form used in generating the grid).
subroutine LFFC(a_st,rho_st,aex_st,g_st_ind,pars,x_guess_gp,bcsf_loss,CTS_ind,x_pol_all_unf,grids_tmp,SS,maxatt_iter)
    Use nag_library, Only: e04jyf, nag_wp, e04jcf, e04jcp,e05jaf,e05jbf,E05JBK,E04ABA

    !Values of state variables as saved in the grid (hence the suffix st).
    real(dp), dimension(I_par), intent(in) :: a_st
    real(dp), dimension(I_par-1), intent(in) :: rho_st
    real(dp), intent(in) :: aex_st
    integer, intent(in) :: g_st_ind
    type(par), intent(in) :: pars

    real(dp), intent(inout), dimension(ncvar_par) :: x_guess_gp !initial guess at the grid point in vector form. On input it contains
    !some naive initialization, which will be overwritten on output by the guess found by subroutine LFFC
    type(grids_type), intent(in) :: grids_tmp !This will be used if we are doing interpolation from previous stage of CTS

    type(ss_type), dimension(0:M_par), intent(in) :: SS !steady-state. The first element (index 0) corresponds to steady state when
    !government expenditure is kept forever at the expected (long-run SS probabilities) government expenditure.
    !Element with index m contains the same, but the expenditure corresponds to shock realization with index m.
    !This is useful for getting initial guess in first period, and for debugging.

    real(dp), intent(out) :: bcsf_loss !minimum loss (that will either be zero if strictly feasible
    !point found or the loss at the minimum loss point that was found).

    integer, intent(in) :: maxatt_iter

    integer, intent(in) :: CTS_ind !stage of CTS algorithm.
    real(dp), dimension(:,:,:,:,:,:), allocatable, intent(in) :: x_pol_all_unf !Policy function from previous stage of CTS algorithm

    real(dp), dimension(3*M_par - 1) :: bcsf_x !We save the best choice found so far here.
    real(dp), dimension(ncvar_par) :: x_interp !initial guess obtained using interpolation.
    integer :: getavfail = 0
    logical :: constr_ok
    real(dp) :: loss

    real(dp), dimension(M_par,I_par) :: l,a_pr
    !real(dp), dimension(1,I_par) :: theta_prime
    real(dp), dimension(M_par,1) :: rho_pr !rho_prime (next-period rho for all states)
    real(dp), dimension(M_par,1) :: aex_pr !mu-adjusted external debt of country 1

    real(dp) :: c_ss_tot !total consumption in steady state
    real(dp) ::  c_tot_perm,rho_perm,l_perm, rnd_num(3) !permutations and 3 random numbers used in generating initial guess

    !Local choice variables: gp stands for grid point.
    real(dp), dimension(M_par,I_par) :: c_gp

    integer :: inform = 0 !local variable for testing loss_func, delete later

    !integer :: att_max !maximum number of attempts (stop when we find a strictly feasible choice)
    !if we do not find a feasible choice use the best of the att_max attempts in terms of loss (not adjusted - see
    !subroutine check_constr for loss adjustment meaning).
    integer :: att_counter !attempt counter (number of attempts tried).

    !real(dp), dimension(ncvar_par) :: xc !choice variables (M_par x (I_par+1) matrix in row form)
    real(dp) :: fc !function value in maximization
    !The following two arrays are arrays of reals and integers which can be used to supply information
    !to the maximized function as an alternative to COMMON block.
    real(dp), dimension(1) :: ruser = 0.0_dp !Contains no information, only defined for purposes of using a NAG subroutine
    integer, dimension(3) :: iuser = 0 !will be used to pass what_binds to some subroutines.
    integer :: ifail = 0
    integer :: nf
    !Bounds:
    real(dp), dimension(ncvar_par) :: bl = 0.00001_dp !lower bound for labour and consumption choices
    real(dp) :: fmax

    integer :: opt_subroutine_local

    !declarations of copies of variables used in the common block
    real(dp), dimension(I_par) :: a_st_cp
    real(dp), dimension(I_par-1) :: rho_st_cp
    real(dp) :: aex_st_cp
    integer :: g_st_ind_cp

    integer, dimension(4) :: j_guess


    real(dp) :: loss_adj_min !minimized adjusted (to account for strictly satisfying some constraints) loss

    !temporary debug variable
    logical :: use_new_guess

    !The following common blocks are for passing data to loss_func (which needs to have a special interface
    !for use with NAG). We create local copies of the relevant variables. The reason for this is that
    !we pass some of these to LFFC so we might need to use the common block one level higher
    !instead of doing that, i.e., in subroutine get_c_guess_img. But this would not be very clean and would possibly
    !yield errors (no intent(in) checks). .
    !Creating a copy of the data is almost costless (small objects and the algorithm takes much more time in the iterative
    !minimization of loss_func.
    common /loss_func_bloc/ a_st_cp,rho_st_cp,aex_st_cp,g_st_ind_cp
    !The parameters are passed to the subroutine using a pointer from a higher level.


    !we do not use j_guess here.
    j_guess = 0


    bcsf_loss = max_loss !initialize loss at the maximum value permitted
    bcsf_x = 0.001_dp !default initial choice.

    !If we have the policy from previous stage of CTS algorithm, use interpolation to obtain initial guess.
    !If it is feasible, we are done. If not, compare the loss of this point with the initial guess
    !obtained manually here, and choose the better of the two (in the sense of lower loss) to continue.
    !If CTS_ind > 1, we should have the policy function from the previous iteration (coarser grid).
    !Use this to obtain initial guess by interpolation.
    if(CTS_ind>1 .and. pars%CTS_C_interp) then
        !Get the initial guess using interpolation, save it in x_interp
        select case(pars%x_interpolation_mode)
            case(1) !Quadrilinear interpolation
        call interp_x(x_pol_all_unf,[a_st(1),a_st(2),rho_st(1),aex_st],j_guess,g_st_ind,x_interp,grids_tmp%a1_gr,&
            grids_tmp%a2_gr,grids_tmp%rho_gr,grids_tmp%aex_gr,size(grids_tmp%a1_gr),size(grids_tmp%rho_gr),size(grids_tmp%aex_gr))
            case(2) !Local Shepard's interpolation
        call interp_x_shep(x_pol_all_unf,[a_st(1),a_st(2),rho_st(1),aex_st],j_guess,g_st_ind,x_interp,grids_tmp%a1_gr,&
            grids_tmp%a2_gr,grids_tmp%rho_gr,grids_tmp%aex_gr,size(grids_tmp%a1_gr),size(grids_tmp%rho_gr),size(grids_tmp%aex_gr),&
            pars%shep_n,pars%shep_norm,pars%shep_norm_val,pars%shep_floor)
        end select

        !Get all the other variables and loss associated with the interpolated choice
        call get_av(pars,x_interp,a_st,rho_st,aex_st,g_st_ind,getavfail,c_gp,l,a_pr,rho_pr,aex_pr,.false.)
        !If this point violates the non-negativity constraints, then it is a disastrous initial guess,
        !so we ignore it and proceed as in the first stage of CTS (manually constructing the initial guess)
        if(getavfail > 0) go to 102

        !Check constraints.
        call check_constr(pars,c_gp,l,a_pr,rho_pr,aex_pr,constr_ok,loss,.false.,.false.)

        bcsf_loss = loss !save best choice so far (need this if we return now, it's part of output)
        !It's the first time that we check it, so it's definitely not worse than the current value
        !which is set at max_loss, so we don't need to check whether it's better than current bcsf_loss).

        !If this is strictly feasible, we are done (save the choice, loss, and return). Loss should be already zero.
        if(constr_ok) then
            !the choice already is in x_guess_gp
            x_guess_gp = x_interp
            return
        else
            bcsf_x = x_interp !save the point as the best choice so far
        end if
102  end if
    !End of the branch of program getting guess by interpolation from last CTS stage policy function
    !If we got here, the initial guess is not feasible. Then we try to construct an initial guess
    !manually as in the first stage of CTS, and try if that is feasible. If neither is feasible,
    !we pick the better of the two in terms of smaller loss, and we initiate improvement from there.

    !x_guess_gp on entry contains a naive initial guess which does not take into account the
    !state variables. In particular, it does not take into account state rho, and this means that
    !if we use steady-state in which country 1 and country 2 consumption in all states but the last
    !are identical, then consumption in country 2 state M does all the adjustment work, and the
    !choice will usually be far from optimum (and due to non-convexities, if we start with
    !such a guess, we might never converge to the optimum).
    !We construct a better guess by taking rho into account, and by adjusting labour supply so that
    !we guarantee that labour supply in country 2 will never be negative.
    !We save this better initial guess into x_guess_init. If it is not feasible (loss > 0), we
    !then try random permutations of this guess.

    !We want to preserve total consumption the same as in autarky steady-state, just not
    !equal division, but such that the ratio of marginal utilities of consumptions is equal to rho_st
    !in every state (so the initial guess will presrve rho in every state, and consumption of country 2
    !in the last state will not fluctuate wildly).

    !(NOTE: unless the grid for rho is very narrow, taking rho_st into account is extremely important,
    !otherwise we often get a very bad initial guess)


    !New construction of initial guess - works for M_par = 2 only. For each shock realization, we have precomputed the total consumption.
    if(M_par == 2 .and. pars%LFFC_use_new_guess) then
    !The old version of the code without adjustment is still used if M_par > 2. I refer to this here.

    !In every state, get the consumption of country 1 and 2 s.t. the ratio of marginal utilities
    !is equal to the ratio of MU last period (the state variable).

    !The consumption in country 1, state 1
    x_guess_gp(3) = SS(1)%c_ss_tot/(pars%mass(1,1) + pars%mass(1,2)*rho_st(1) ** (1.0_dp/pars%sigma_par))
    !The consumption in country 1, state 2
    x_guess_gp(4) = SS(2)%c_ss_tot/(pars%mass(1,1) + pars%mass(1,2)*rho_st(1) ** (1.0_dp/pars%sigma_par))

    !Consumption of country 2 in state 1
    x_guess_gp(5) = rho_st(1) ** (1.0_dp/pars%sigma_par) * x_guess_gp(3)



    !Now adjust these by the pre-computed adjustment (the total consumption in state 1 will be
    !preserved, but the ratio will be wrong. Also, put a floor on consumption of each agent
    !as 0.01. If in doubt about what these mean, check the Onenote note called 'Improving
    !initial guess'
    x_guess_gp(3) = max(x_guess_gp(3) + SS(1)%kappa_g2ming1,0.01_dp)
    x_guess_gp(4) = max(x_guess_gp(4) + SS(2)%kappa_g2ming1,0.01_dp)
    x_guess_gp(5) = max(x_guess_gp(5) - SS(1)%kappa_g2ming1,0.01_dp)


    !Get the labour supply of country 1. For now assume that the division of production is the same (perhaps
    !we could do slightly better than this but it is not obvious how).
    !state 1 (consumption is x_guess_gp(3)
    x_guess_gp(1) = (pars%mass(1,1)*(x_guess_gp(3) + pars%G(1,1)) + &
              pars%mass(1,2)*(x_guess_gp(5) + pars%G(1,2) ))/(&
              pars%mass(1,1)*pars%theta_0(1,1) + pars%mass(1,2)*pars%theta_0(1,2) )

    !In state 2, because we use the rho constraint to obtain c(2,2), the actual
    !aggregate private consumption will be different than the one in SS which was
    !pre-computed. However, the difference should not be great - so for the purposes
    !of obtaining an initial guess, just take the private consumption from steady state
    !if shock realization index is forever 2. We could also use SS(1) instead of x_guess_gp
    !when computing the labour supply in state 2.
    x_guess_gp(2) = (SS(2)%c_ss_tot + pars%mass(1,1)*pars%G(2,1) + pars%mass(1,2)* pars%G(2,2))/(&
              pars%mass(1,1)*pars%theta_0(1,1) + pars%mass(1,2)*pars%theta_0(1,2) )

    else !For general M_par, I do not currently adjust the initial guess to take into account
    !government expenditures in a more sophisticated manner...
    c_ss_tot = pars%mass(1,1)*x_guess_gp(M_par + 1) + pars%mass(1,2)*x_guess_gp(2*M_par + 1)
    !This assumes functional form of utility function. For consistency, it should be defined in mod_utility as
    !other expressions which depend on functional form.
    !The following elements contain consumption of country 1 in all states
    x_guess_gp(1+M_par:2*M_par) = c_ss_tot/(pars%mass(1,1) + pars%mass(1,2)*rho_st(1) ** (1.0_dp/pars%sigma_par) )
    !Now copy labour supply of country 2 (only if M_par > 1, otherwise we never save that)
    if(M_par > 1) then
    x_guess_gp(1+2*M_par:3*M_par-1) = rho_st(1) ** (1.0_dp/pars%sigma_par) * x_guess_gp(1+M_par)
    end if
    !Finally get labour supply. Choose it in such a way that labour supply (actual l, not production in a country pi * l) will
    !be equal in both countries.
    x_guess_gp(1:M_par) = reshape( (pars%mass(1,1)*(x_guess_gp(1+M_par) * ones_col + pars%G(:,1:1) ) + &
              pars%mass(1,2)*(x_guess_gp(1+2*M_par) * ones_col + pars%G(:,2:2) ))/(&
              pars%mass(1,1)*pars%theta_0(1,1) + pars%mass(1,2)*pars%theta_0(1,2) ) ,[M_par])
    end if


    !First, check whether the initial guess x_guess_0 is feasible (and compute the associated value of loss function).
    !If it is feasible, we are done (feasibility is the only requirement). The subroutine returns the current guess x_guess_init.

    !Call subroutine get_av to get all missing variables and next-period state for the choice variables.
    call get_av(pars,x_guess_gp,a_st,rho_st,aex_st,g_st_ind,getavfail,c_gp,l,a_pr,rho_pr,aex_pr,.false.)

    !If non-negativity of consumption or labour is violated, then current return is not defined. Also, it is probably pointless
    !to try any local improvements. Just go to the next part of the algorithm which initializes improvements starting at random
    !initial guesses.
    if(getavfail > 0) then
    !Most likely this was due to the intial guess adjustment for g. Just use unadjusted version in this case.
    !(this will not be used as long as any of the permutations is better, but we need to initialize this, otherwise
    !we might get crashes/errors elsewhere).
    !The consumption in country 1, state 1
    x_guess_gp(3) = SS(1)%c_ss_tot/(pars%mass(1,1) + pars%mass(1,2)*rho_st(1) ** (1.0_dp/pars%sigma_par))
    !The consumption in country 1, state 2
    x_guess_gp(4) = SS(2)%c_ss_tot/(pars%mass(1,1) + pars%mass(1,2)*rho_st(1) ** (1.0_dp/pars%sigma_par))
    !Consumption of country 2 in state 1
    x_guess_gp(5) = rho_st(1) ** (1.0_dp/pars%sigma_par) * x_guess_gp(3)

    !The labour supplies which we got should already be consistent with this solution (i.e., we
    !should never get a negative labour supply) if we are using the new algorithm. If we are using
    !the old algorithm, at this stage this error should never happen.

        goto 101
    end if

    !Compute the unadjusted loss (so the last argument of check_constr is false). This is used for ranking
    !the points because feasibility has priority. In any case, whether we should the adjusted or unadjusted
    !loss shouldn't make much difference (they are closely related and in the case where ranking reversed,
    !we would still expect the points to be close and the maximization in later stages of the algorithm should
    !thus converge to a similar point even in the presence of nonconvex feasible sets.).
    call check_constr(pars,c_gp,l,a_pr,rho_pr,aex_pr,constr_ok,loss,.false.,.false.)

    !This is the first point that we tried. Therefore, best choice so far loss is the current loss, and the best choice
    !so far guess (x) is the current guess - that was already initialized at the start. This may change if
    !we start using the interpolation.
    !See if this is better than the bcsf (this could either be a maximum number, or the loss associated with interpolation
    !If this is the best point so far, then set bcsf_x to x_guess_gp, and update bcsf_loss. Otherwise bcsf_x remains
    !x_interp
    if(loss <= bcsf_loss) then
        bcsf_x = x_guess_gp !and initialize best choice so far
        bcsf_loss = loss
    end if

    !(fc_found in place of rofv because if we haven't found a feasible choice yet, we are potentially interested in
    !non-feasible ones too - then we pick the one with the minimum loss - so we need to check all constraints and compute
    !the loss function)

    !If this is strictly feasible, we are done (save the choice, loss, and return).
    if(constr_ok) then
        !the choice already is in x_guess_gp so we don't need to copy it.
        return
    end if

    !If the choice is not feasible, then proceed with improving the initial
    !guess using an optimization approach (loss minimization), starting with random initial guesses
    !In the first attempt, start directly at the choice, then proceed with randomization.

!Create copies of states for use in common bloc (to be passed to loss_func)
!(Do this before starting the cycle)
    !Before getting an initial guess of consumption and proceeding with the loss minimization, create
    !a copy of all the variables in the common blocks shared by subroutine loss_func.
101 a_st_cp = a_st
    rho_st_cp = rho_st
    aex_st_cp = aex_st
    g_st_ind_cp = g_st_ind

    att_counter = 0 !initialize number of attempts
    iuser = 0 !initialize (this variable needs to be defined for the NAG routines but is not used)

    !Generate a random permutation of the initial guess x_init, and save it in x_guess%x (as starting point
    !for improving algorithm). This is just like generating the initial guess x_init (with the added permutations)
    !and it could be at some point put in a subroutine but it's low priority.
    !If this is not the first attempt (we initially want to try to improve the naive initial guess), permutate
150 if(att_counter > 0) then
        call random_number(rnd_num)

        if(M_par /= 2 .or.  (.not. pars%LFFC_use_new_guess)) then
        !use the old algorithm for permutation (which works for arbitrary M)

        !Scale these numbers so they are in proper intervals (first in (1-LFFC_rnd_a,1+LFFC_rnd_a), etc.)
        !Apply multiplicative shocks
        c_tot_perm = c_ss_tot * (1.0_dp - pars%LFFC_rnd_a + rnd_num(1)*2.0_dp*pars%LFFC_rnd_a)
        rho_perm = rho_st(1) * (1.0_dp - pars%LFFC_rnd_b + rnd_num(2)*2.0_dp*pars%LFFC_rnd_b)
        l_perm = 1.0_dp * (1.0_dp - pars%LFFC_rnd_c + rnd_num(3)*2.0_dp*pars%LFFC_rnd_c)

        !Now compute the initial guess and save it in x_guess_gp. Similarly as we compute x_guess_init.
        !The only difference is that we use the c_tot_perm instead of the steady-state consumption,
        !rho_perm instead of the state rho_st, and l_perm (ratio of labour supplies) which deviates from
        !1 (not from steady-state ratio - that might not be closely related to the consumptions anyway
        !once rho is taken into account, which is not taken into account when computing steady-state labour supply).

        x_guess_gp(1+M_par:2*M_par) = c_tot_perm/(pars%mass(1,1) + pars%mass(1,2)*rho_perm ** (1.0_dp/pars%sigma_par) )
        !Now copy labour supply of country 2 (only if M_par > 1, otherwise we never save that)
        if(M_par > 1) then
            x_guess_gp(1+2*M_par:3*M_par-1) = rho_perm ** (1.0_dp/pars%sigma_par) * x_guess_gp(1+M_par)
        end if
        !Get labout supply. Choose it in such a way that labour supply ratio (actual l, not production in a country pi * l) will
        !be l_perm (l2/l1 = l_perm) as opposed to 1 which we used when constructing the initial guess (so l_perm > 1 means
        !that household in country 2 work more).
        x_guess_gp(1:M_par) = reshape( (pars%mass(1,1)*(x_guess_gp(1+M_par) * ones_col + pars%G(:,1:1) ) + &
            pars%mass(1,2)*(x_guess_gp(1+2*M_par) * ones_col + pars%G(:,2:2) ))/(&
            pars%mass(1,1)*pars%theta_0(1,1) + pars%mass(1,2)*pars%theta_0(1,2)*l_perm ) ,[M_par])

        else !use the new algorithm

        !warning - the meanings of c_tot_perm and rho_perm are different (they are now distributed around 1 symmetrically,
        !rather than being already the permutated variables)! [different compared to he previous algorithm.]
        !Scale these numbers so they are in proper intervals (first in (1-LFFC_rnd_a,1+LFFC_rnd_a), etc.)
        !Apply multiplicative shocks
        c_tot_perm = (1.0_dp - pars%LFFC_rnd_a + rnd_num(1)*2.0_dp*pars%LFFC_rnd_a)
        !This increases/decreases the total consumption in all states

        !There are two options. The default one (for positive values of LFFC_rnd_b is that
        !kappa is subjected to a multiplicative shock, then is used as before - this leads
        !to sometimes the adjustment being smaller and sometimes greater but always around the
        !initial guess of adjustment). If LFFC_rnd_b is negative, then the multiplicative
        !shock to kappa is drawn from (0,1) interval. The interpretation is that kappa
        !in parameter file provides an upper bound on adjustment, and the adjustment
        !is varied between 0 and this upper bounds.

        if(pars%LFFC_rnd_b < 0.0_dp) then
            rho_perm = rnd_num(2)
        else
            rho_perm = (1.0_dp - pars%LFFC_rnd_b + rnd_num(2)*2.0_dp*pars%LFFC_rnd_b)
        end if

        !I still call this rho_perm so I don't have to use another variable but it is something
        !different (it does not contain rho, it's just number uniformly drawn from the interval
        !(1-LFFC_rnd_b,1+LFFC_rnd_b).
        !rho_perm should perhaps be called something like kappa_g2ming1_perm, but it is the same
        !in the sense that it permutates the deviation from rho_st (in terms of marginal utilities)
        l_perm = (1.0_dp - pars%LFFC_rnd_c + rnd_num(3)*2.0_dp*pars%LFFC_rnd_c)

        !c_tot_perm permutates the levels of total consumption in all states (actually in the M=2
        !case it will be in state 1 only exactly. The state 2 consumption will not be quite the same.

        !The consumption in country 1, state 1
        x_guess_gp(3) = SS(1)%c_ss_tot*c_tot_perm/(pars%mass(1,1) + pars%mass(1,2)*rho_st(1) ** (1.0_dp/pars%sigma_par))
        !The consumption in country 1, state 2
        x_guess_gp(4) = SS(2)%c_ss_tot*c_tot_perm/(pars%mass(1,1) + pars%mass(1,2)*rho_st(1) ** (1.0_dp/pars%sigma_par))

        !Consumption of country 2 in state 1 (start off with the assumption that the rho constraint is satisfied state by state)
        x_guess_gp(5) = rho_st(1) ** (1.0_dp/pars%sigma_par) * x_guess_gp(3)


        !Now adjust these by the pre-computed adjustment (the total consumption in state 1 will be
        !preserved, but the ratio will be wrong. Also, put a floor on consumption of each agent
        !as 0.01. If in doubt about what these mean, check the Onenote note called 'Improving
        !initial guess'
        x_guess_gp(3) = max(x_guess_gp(3) + SS(1)%kappa_g2ming1*rho_perm,0.01_dp)
        x_guess_gp(4) = max(x_guess_gp(4) + SS(2)%kappa_g2ming1*rho_perm,0.01_dp)
        x_guess_gp(5) = max(x_guess_gp(5) - SS(1)%kappa_g2ming1*rho_perm,0.01_dp)


        !Get the labour supply of country 1. For now assume that the division of production is the same (perhaps
        !we could do slightly better than this but it is not obvious how).
        !state 1 (consumption is x_guess_gp(3)
        x_guess_gp(1) = ((pars%mass(1,1)*(x_guess_gp(3) + pars%G(1,1)) + &
                pars%mass(1,2)*(x_guess_gp(5) + pars%G(1,2) ))/(&
                pars%mass(1,1)*pars%theta_0(1,1) + pars%mass(1,2)*pars%theta_0(1,2) ))*l_perm

        !In state 2, because we use the rho constraint to obtain c(2,2), the actual
        !aggregate private consumption will be different than the one in SS which was
        !pre-computed. However, the difference should not be great - so for the purposes
        !of obtaining an initial guess, just take the private consumption from steady state
        !if shock realization index is forever 2. We could also use SS(1) instead of x_guess_gp
        !when computing the labour supply in state 2.
        x_guess_gp(2) = ((SS(2)%c_ss_tot + pars%mass(1,1)*pars%G(2,1) + pars%mass(1,2)* pars%G(2,2))/(&
                pars%mass(1,1)*pars%theta_0(1,1) + pars%mass(1,2)*pars%theta_0(1,2) ))*l_perm


        !NOTE: IT can happen that the labour supply of country 2 (after permutation is actually negative.
        !This is because of the problem described in the previous paragraph. It should almost never happen and if it does
        !it will just lead to getavfail at some point, and a return capped at maximum loss. So such points should
        !never be chosen. IF they are, then we can fix this.
        !
        !We could fix this by solving for the missing element of the consumption matrix using the
        !rho constraint. Then we could get the aggregate consumption in state 2, and we could check
        !that for the chosen labour supply in country 1 (after permutation), the labour supply in country 2 is positive.
        !But it is quite cumbersome and it should not matter, because due to the low
        !returns, such points should never be chosen in any region of state space (we just need to
        !increase LFFC_attmax if this happens).
        !
        !The easiest thing to do is to check whether this happens and if yes, then
        !take a guess in which we permutate only the total consumption level (we need to recompute labour supply
        !too).
        end if !end of using the branch for using the old algorithm.

        !Now we got a random guess (wherever it cam from) - we need to check that it is ok.
        !IF not - then we need to fix the guess.

        !Check the initial guess:
        call get_av(pars,x_guess_gp,a_st,rho_st,aex_st,g_st_ind,getavfail,c_gp,l,a_pr,rho_pr,aex_pr,.true.)

        !The same as above but we skip the step where we permutate the rhos and labour supply! This should
        !guarantee that get_av_fail doesn't happen. Actually already skipping the rho permutation should
        !be enough, but it doesn't matter. This should happen in very rare cases.

        !OK, there is a big problem here - the following should not be done always (that just wastes
        !the permutation because we are effectively permutating in one dimension only. BUT now that I put it
        !into the following if condition, the program crashes (negative consumption sometimes apparently).
        !I need to get to the bottom of this.

        if(getavfail > 0) then

        att_counter = att_counter + 1
        if(att_counter < maxatt_iter) goto 150

        !In the last attempt if getavfail we fix the guess. This should never be necessary because we should
        !already have a guess which works from this respect (non-neg constraint) but it is convenient
        !(otherwise we would have to handle some errors a bit more carefully down the line) and it is
        !not very costly.
        x_guess_gp(3) = SS(1)%c_ss_tot*c_tot_perm/(pars%mass(1,1) + pars%mass(1,2)*rho_st(1) ** (1.0_dp/pars%sigma_par))
        x_guess_gp(4) = SS(2)%c_ss_tot*c_tot_perm/(pars%mass(1,1) + pars%mass(1,2)*rho_st(1) ** (1.0_dp/pars%sigma_par))
        !Consumption of country 2 in state 1 (start off with the assumption that the rho constraint is satisfied state by state)
        x_guess_gp(5) = rho_st(1) ** (1.0_dp/pars%sigma_par) * x_guess_gp(3)

        x_guess_gp(1) = ((SS(1)%c_ss_tot*c_tot_perm + pars%mass(1,1)*(pars%G(1,1)) + &
                pars%mass(1,2)*( pars%G(1,2) ))/(&
                pars%mass(1,1)*pars%theta_0(1,1) + pars%mass(1,2)*pars%theta_0(1,2) ))
        x_guess_gp(2) = ((SS(2)%c_ss_tot*c_tot_perm + pars%mass(1,1)*pars%G(2,1) + pars%mass(1,2)* pars%G(2,2))/(&
                pars%mass(1,1)*pars%theta_0(1,1) + pars%mass(1,2)*pars%theta_0(1,2) ))
        end if

    else !If this is the first attempt
        x_guess_gp = bcsf_x !this is either the initially constructed guess (unpermutated), or the interpolated
        !policy value from previous stage of CTS.
    end if

    att_counter = att_counter + 1 !increase the counter of attempts made

    !In the maximization - be very careful not to evaluate current return when getavfail > 0!!!
    !(instead give a pre-defined value (max_loss in hard parameter), and ideally return control
    !(check documentation of NAG subroutine) so we don't waste time trying to locally improve
    !a hopeless solution. Instead start with a new one.



    !The initial guess is constructed in such way that negative consumption or labour supply should never happen)
    !(do not check at this stage to save time - but we could use get_av with the last parameters false).

    ifail = 1 !silent exit (we do not want the program to stop when soft errors happen - for example sometimes
    !it is bound to happen that e04jcf does not find the optimum within the limit of iterations).

    opt_subroutine_local = 1 !as opposed to pars%opt_subroutine
    !(1 works fine, but we sometimes want to use 2, etc. in simulations). Rerwiting this choice here
    !avoids crashes due to a mistake in the parameter file.
    select case(opt_subroutine_local)
        case(1) !Quadratic approximation method
            call e04jcf(loss_func,ncvar_par,ncvar_par*2+1,x_guess_gp,bl,pars%x_bu,0.05_dp,&
                pars%rhoend,e04jcp,pars%LFFC_maxiter,loss_adj_min,nf,iuser,ruser,ifail)

        case(2) !Newton-like method (leads to crashes on HPC) - rather use the quadratic approx. method
        !(maybe try to use this one at some point later).
            !call e04jyf(M_par*I_par-1,0,CPC_ret,bl,bu,xc,fmax,iw,liw,w,lw,iuser,ruser,ifail)
            if(this_image()==1) then
                write(*,*) 'Error: need to generalize LFFC for optim_subroutine = 2'
                error stop
            end if
        case default
            write(*,*) 'Error in subroutine LFFC - wrong value of optim_subroutine'
            error stop
    end select

    !Now compute the loss of the point (unadjusted) to evaluate how good the point is in terms of feasibility
    !because after all we are firstly looking for strictly feasible points), and being strictly inside bounds
    !(which the adjusted loss rewards) is just a bonus. This is the same as what we do a bit earlier in the subroutine.
    !(explanation there)

    call get_av(pars,x_guess_gp,a_st,rho_st,aex_st,g_st_ind,getavfail,c_gp,l,a_pr,rho_pr,aex_pr,.false.)
    if(getavfail > 0 .and. att_counter < maxatt_iter) then
        !This should perhaps be handled better - the fact that there is no crash
        !coming from here depends on the assumption that there are enough attempts,
        !so that there will (at every grid point) always be a point which satisfies
        !at least the non-neg. constraints. This might not be always the case if
        !we are very unlucky (huge number of grid points and low number of attempts).
        att_counter = att_counter + 1 !I am not 100% sure if this is needed but in the worst
        !case we make 1 less attempt (if Ididn't put it here maybe there could be an infinite
        !cycle sometimes?)
        go to 150 !try again with a different initial guess
    end if
    !Compute the loss
    !(if getavfail > 0 set it to maxloss instead) - this was not handled properly earlier as there was a chance
    !that the program might crash due to passing nagative c or l to check_constr
    if(getavfail > 0) then
        loss = max_loss
    else
        call check_constr(pars,c_gp,l,a_pr,rho_pr,aex_pr,constr_ok,loss,.false.,.false.)
    end if

    !If loss < 0.001_dp,treat the point as feasible, save the loss, and return.
    !(This is only a very slight violation of a constraint, in practice no problem).
    !(above, we were a bit stricter and checked strict feasibility, i.e., loss exactly zero). But
    !if nothing worked so far, we should be ready to compromise and accept a slightly infeasible
    !guess here (loss = 0.001 corresponds, given the strict loss function, to the slightest violation).
    if(loss < 0.001_dp) then
        bcsf_loss = loss !save the loss of best choice so far and return (at this point x_guess_gp%x contains
        !the choice for which this minimal loss was attained)
        return
    else
        !If the point is better than all others so far, save it and its loss
        if(loss < bcsf_loss) then
            bcsf_loss = loss
            bcsf_x = x_guess_gp
        end if
        !If we have not reached the maximum number of attempts, try again to find a better point
        if(att_counter < maxatt_iter) go to 150
    end if

!If no feasible choice found (so no return call yet), save the least infeasible choice and return loss so we can compute average loss
!and get an idea of how bad the violation of constraints is
x_guess_gp = bcsf_x
!(the returned loss is bcsf_loss so there is no need to assign that to anything)

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


!The subroutine checks whether file stop_now.txt exists in the folder results/folder_name.
!If yes this is then used to stop execution. This is useful so that we can send a signal to
!stop VFI (or simulation) and take what's been done so far (rather than wait for
!runtime limit to run out). Especially in VFI sometimes the program might converge faster
!than expected, yet thr stopping rule might never be quite satisfies.
subroutine  stop_now_check(folder_name,stop_now)
    character*(*), intent(in) :: folder_name !folder_name is the output folder (in which log is saved)
    logical, intent(out) :: stop_now
    character*256 :: shell_cmd,file_path

    logical :: file_exists

    inquire(file=('results/'//trim(folder_name)//'/stop_now.txt'), exist=file_exists)
    if(file_exists) then
        stop_now = .true.
    else
        stop_now = .false.
    end if

end subroutine stop_now_check

!subroutine load_N loads the number of gridpoints. This is then used to allocate arrays for
!loading the other variables. M_par_loc is loaded only to check that the value
!function in the file was generated for the same value of M_par, otherwise it makes no sense to use
!the value function loaded from file.
subroutine load_N(N_a,N_rho,N_aex,input_folder,folder_name,M_grid)
    integer, intent(out) :: N_a,N_rho,N_aex

    character*(*), intent(in) :: input_folder,folder_name !folder_name is the output folder (in which log is saved)
    character*256 :: shell_cmd,file_path

    integer, intent(in) :: M_grid !This tells the program what the M_grid is (to check that this
    !is the same as M_par_loc) when loading this.

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

    inquire(file=('results/'//trim(input_folder)//'/N_aex.out'), exist=file_exists)
    if (.not. file_exists) then
        call double_output(folder_name,'Error: File N_aex.out in input_folder not found.')
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

    file_path = 'results/'//trim(input_folder)//'/N_aex.out'
    open(unit=20, file = file_path, status = 'old')
    read(unit=20, fmt=*)  N_aex
    close(unit = 20)

    file_path = 'results/'//trim(input_folder)//'/M_par_loc.out'
    open(unit=20, file = file_path, status = 'old')
    read(unit=20, fmt=*)  M_par_loc
    close(unit = 20)

    if(M_par_loc /= M_grid) then
        write(*,*) 'Error: M_par_loc does not equal M_grid.'
        write(*,*) 'Value function loaded from file was computed for different M_par. It is useless as initial guess.'
        write(*,*) 'Possible cause is different setting of parameter IID_shocks for the two calibrations.'
        error stop
    end if

end subroutine load_N

!subroutine load_V loads the value function and grids. Interpolation is used. Due to numerical errors,
!it is easiest to use interpolation every time even if the grid were identical.
subroutine load_V(V_vect,a1_gr,a2_gr,rho_gr,aex_gr,input_folder,folder_name)
    real(dp), dimension(:), intent(out) :: V_vect
    real(dp), dimension(:), intent(out) :: a1_gr,a2_gr,rho_gr,aex_gr

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
    inquire(file=('results/'//trim(input_folder)//'/aex_gr.out'), exist=file_exists)
    if (.not. file_exists) then
        call double_output(folder_name,'Error: File aex_gr.out in input_folder not found.')
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

    file_path = 'results/'//trim(input_folder)//'/aex_gr.out'
    open(unit=20, file = file_path, status = 'old')
    read(unit=20, fmt=*)  aex_gr
    close(unit = 20)

end subroutine load_V


subroutine load_policy(x_pol_all,input_folder,folder_name)
    real(dp), dimension(:,:), intent(out) :: x_pol_all

    character*(*), intent(in) :: input_folder,folder_name !folder_name is the output folder (in which log is saved)
    character*256 :: file_path

    !no checking whether the files exist (it was done before calling this subroutine).

    !Reading all the variables
    file_path = 'results/'//trim(input_folder)//'/x_pol_all.out'
    open(unit=20, file = file_path, status = 'old', form='binary')
    read(unit=20)  x_pol_all
    close(unit = 20)

end subroutine load_policy


!The following subroutine computes maximum feasible consumption going to a country
!for given productivity and government expenditure in levels.
subroutine max_feas_cons(l_max,mass,theta_lvl,g_lvl,max_cons)
    real(dp), dimension(1,I_par), intent(in) :: l_max,theta_lvl,mass
    real(dp), dimension(M_par,I_par) ,intent(in) :: g_lvl

    real(dp),dimension(M_par,I_par), intent(out) :: max_cons

    integer :: m_ind

    !Reminder: This is maximum consumption by a country
    do m_ind = 1,M_par
        max_cons(m_ind,:) = sum(mass*theta_lvl*l_max) - mass(1,1)*g_lvl(m_ind,1) - mass(1,2)*g_lvl(m_ind,2)
    end do

    !Now divide it by the size of country to get per-capita maximum consumption (that's the actual choice variable)
    max_cons(:,1) = max_cons(:,1)/mass(1,1)
    max_cons(:,2) = max_cons(:,2)/mass(1,2)

    if(I_par > 2) then
        if(this_image()==1) write(*,*) 'Error: subroutine max_feas_cons needs to be generalized for I_par > 2.'
        error stop
    end if
    sync all

end subroutine max_feas_cons


!Subroutine save_V saves the value function and the grid for which it was generated.
!It also saves the number of grid points. This is useful when we load the
!value function from file. M_par is saved just for checking purposes - in practice
!we will never want to use a value function computed for different M.
subroutine save_V(V_vect,grids,N_a,N_rho,N_aex,M_par_loc,folder_name)
    real(dp), dimension(:), intent(in) :: V_vect
    type(grids_type), intent(in) :: grids
    integer, intent(in) :: N_a,N_rho,N_aex,M_par_loc !M_par_loc is just a local copy

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

    file_path = 'results/'//trim(folder_name)//'/aex_gr.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  grids%aex_gr
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/N_a.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  N_a
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/N_rho.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  N_rho
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/N_aex.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  N_aex
    close(unit = 20)

    !This might be less than M_par (for example if IID_shocks used,
    !this will be 1 irrespective of valur of M_par).
    file_path = 'results/'//trim(folder_name)//'/M_par_loc.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  M_par_loc
    close(unit = 20)
end subroutine save_V

!Subroutine save_par saves some parameters in the parameter file at runtime in folder
!results/folder_name/pars. This could be useful in the future if we want to
!compare results obtained under different parameters in matlab (and loading
!these directly from the parameter file which contains the same information would be
!cumbersome). Only the parameters related to the environment are saved, parameters
!related to technical stuff (such as number of iterations, etc. are not - these can
!still be recovered from the text file and there is no use for them in Matlab in
!generating plots).
subroutine save_par(pars,folder_name)
    type(par), intent(in) :: pars

    character*(*), intent(in) :: folder_name
    character*256 :: shell_cmd,file_path

    !Saving all the parameters

    file_path = 'results/'//trim(folder_name)//'/pars/T_sim.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  pars%T_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/pars/N_sim.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  pars%N_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/pars/mass.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  pars%mass
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/pars/l_max.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  pars%l_max
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/pars/beta.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  pars%beta
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/pars/theta_0.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  pars%theta_0
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/pars/IID_shocks.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  pars%IID_shocks
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/pars/P.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  pars%P
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/pars/G.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  pars%G
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/pars/s_init_index.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  pars%s_init_index
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/pars/b_init.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  pars%b_init
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/pars/b_ex_init.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  pars%b_ex_init
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/pars/alpha.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  pars%alpha
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/pars/gamma_par.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  pars%gamma_par
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/pars/A_par.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  pars%A_par
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/pars/B_par.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  pars%B_par
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/pars/sigma_par.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*)  pars%sigma_par
    close(unit = 20)

end subroutine save_par

!Subroutine save_policy saves the policy function
subroutine save_policy(x_pol_all,folder_name)
    real(dp), dimension(:,:), intent(in) :: x_pol_all

    character*(*), intent(in) :: folder_name
    character*256 :: shell_cmd,file_path,tmp_string

    file_path = 'results/'//trim(folder_name)//'/x_pol_all.out'
    open(unit=20, file = file_path, status = 'replace', form='binary')
    write(unit=20)  x_pol_all
    close(unit = 20)

end subroutine save_policy

!SWF(c,l,s_current) is expected social welfare, given a matrix of state-contingent
!consumption and labour.
!Each row of c and l corresponds to a state, each column to an agent. s_ind_st is the
!index of last-period shock realization (state). mass and alpha are a vector of
!mass of agents and Pareto weights in the social welfare function.
function SWF(c,l,s_ind_st,P,mass,alpha,A_par,B_par,sigma_par,Gamma_par)
    real(dp) :: SWF

    real(dp), dimension(M_par,I_par), intent(in) :: c,l
    real(dp), dimension(M_par,M_par), intent(in) :: P
    integer, intent(in) :: s_ind_st
    real(dp), dimension(1,2), intent(in) :: mass,alpha
    real(dp), intent(in) :: A_par,B_par,sigma_par,Gamma_par !parameters of the utility function

    real(dp), dimension(M_par,I_par) :: u !Utility of all agents in all states.

    call util(c,l,U,A_par,B_par,sigma_par,Gamma_par) !need to use mod_util for this

    SWF = maxval(matmul(matmul(P(s_ind_st:s_ind_st,1:M_par),U),transpose(mass*alpha)))
    !The maxval is there to convert array of dimension (1,1) to scalar. Perhaps this is not the
    !most efficient way to do this but I don't know how else to do it in a single statement,
    !i.e., without saving the result as (1,1) array and then accesing this array.
end function SWF


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

    integer :: a1_ind,a2_ind,rho_ind,aex_ind,g_ind

    integer :: M_grid

    integer, dimension(M_par,4) :: j_guess

    !Perform checks that the grids and the value functions are compatible in the sense of dimensionality.

    shape_old = shape(V_old)
    shape_new = shape(V_new)

    shape_grid_old(1) = size(grids_old%a1_gr)
    shape_grid_old(2) = size(grids_old%a2_gr)
    shape_grid_old(3) = size(grids_old%rho_gr)
    shape_grid_old(4) = size(grids_old%aex_gr)
    !(grids for shock realization not checked). They should always
    !be the same and we only check the whether the value functions
    !have the same M_par. Also, no interpolation happens over shock
    !realizations - we just assume that shock space is the same here.

    shape_grid_new(1) = size(grids_new%a1_gr)
    shape_grid_new(2) = size(grids_new%a2_gr)
    shape_grid_new(3) = size(grids_new%rho_gr)
    shape_grid_new(4) = size(grids_new%aex_gr)

    if(any(shape_old(1:4)/=shape_grid_old(1:4))) then
        call double_output(folder_name, 'Error: grids_old in subroutine V_from_V incompatible with V_old.')
        error stop
    end if
    if(any(shape_new(1:4)/=shape_grid_new(1:4))) then
        call double_output(folder_name, 'Error: grids_new in subroutine V_from_V incompatible with V_new.')
        error stop
    end if

    if(shape_old(5) /= shape_new(5)) then
        call double_output(folder_name, 'Error: Mismatch of number of gridpoints in V_from_V.')
        error stop
    end if
    M_grid = shape_old(5)
    !This may be different from M_par because of IID_shocks possibility


    !ALSO check whether the ranges of the new grid exceed those of the old grid. If yes,
    !extrapolation is used. For small differences this is fine but a warning should be displayed.
    !(do not assume monotonicy in grids)
    if((minval(grids_new%a1_gr)<minval(grids_old%a1_gr)) .or. &
    (minval(grids_new%a2_gr)<minval(grids_old%a2_gr)) .or. &
    (minval(grids_new%rho_gr)<minval(grids_old%rho_gr)) .or. &
    (minval(grids_new%aex_gr)<minval(grids_old%aex_gr)) .or. &
    (maxval(grids_new%a1_gr)>maxval(grids_old%a1_gr)) .or. &
    (maxval(grids_new%a2_gr)>maxval(grids_old%a2_gr)) .or. &
    (maxval(grids_new%rho_gr)>maxval(grids_old%rho_gr)) .or. &
    (maxval(grids_new%aex_gr)>maxval(grids_old%aex_gr))) then
        call double_output(folder_name,'Warning: in V_from_V, bounds of grid_new exceed bounds of grid_old.')
        call double_output(folder_name,'Extrapolation will be used, dangerous far from bounds.')

        !Also write more details about this in terminal

    write(*,*) '==========Detailed debug informatio==================='
    write(*,*) 'Details of the warning:'
    write(*,*) 'old grids boundaries: '
    write(*,*) 'a1: ',minval(grids_old%a1_gr),maxval(grids_old%a1_gr)
    write(*,*) 'a2: ',minval(grids_old%a2_gr),maxval(grids_old%a2_gr)
    write(*,*) 'rho: ',minval(grids_old%rho_gr),maxval(grids_old%rho_gr)
    write(*,*) 'aex: ',minval(grids_old%aex_gr),maxval(grids_old%aex_gr)
    write(*,*) 'new grids boundaries: '
    write(*,*) 'a1: ',minval(grids_new%a1_gr),maxval(grids_new%a1_gr)
    write(*,*) 'a2: ',minval(grids_new%a2_gr),maxval(grids_new%a2_gr)
    write(*,*) 'rho: ',minval(grids_new%rho_gr),maxval(grids_new%rho_gr)
    write(*,*) 'aex: ',minval(grids_new%aex_gr),maxval(grids_new%aex_gr)


    write(*,*) '======================================================'
    end if

    !In this case we do not have a good initial guess to use (but if we do not initialize this, the program
    !will crash).
    j_guess = 0

    do a1_ind = 1,shape_grid_new(1)
        do a2_ind = 1,shape_grid_new(2)
            do rho_ind = 1,shape_grid_new(3)
                do aex_ind = 1,shape_grid_new(4)
                    do g_ind = 1,M_grid !if shocks are IID, we do not cycle over all grid points.
                    !get the interpolated value and save it in the appropriate element of V_new
                    call interp_V(V_old,[grids_new%a1_gr(a1_ind),grids_new%a2_gr(a2_ind),grids_new%rho_gr(rho_ind),&
                        grids_new%aex_gr(aex_ind)],j_guess(g_ind,:),g_ind,V_new(a1_ind,a2_ind,rho_ind,aex_ind,g_ind)&
                        ,grids_old%a1_gr,grids_old%a2_gr,grids_old%rho_gr,grids_old%aex_gr,&
                        shape_grid_old(1),shape_grid_old(3),shape_grid_old(4))
                end do
            end do
        end do
    end do
end do

write(*,*) 'WARNING: subroutine V_from_V assumes quadrilinear interpolation. We might want to generalize that!!!'


end subroutine V_from_V

!Subroutine VFI_solve_gp solves the value function iteration maximization step at a single
!grid point. This subroutine needs NAG libraries.
subroutine VFI_solve_gp(ind_unfold,grids,pars,N_a,N_rho,N_aex,x_guess_gp,V_new_gp,skip_max,loss_gp,x_AARD_gp,x_MARD_gp,&
    x_AAD_gp,maxatt,gp_ind,runtime_limit_unix,PSO_this_iter)
    Use nag_library, Only: e04jyf, nag_wp, e04jcf, e04jcp,e05jaf,e05jbf,E05JBK,E04ABA,e05saf,E05SXM
    !USE nag_library

    integer, dimension(5), intent(in) :: ind_unfold !vector of indices (to get current state)
    type(grids_type), intent(in) :: grids
    type(par), intent(in) :: pars
    integer, intent(in) :: N_a,N_rho,N_aex !number of gridpoints (could just recover them from the grids).
    real(dp), dimension(ncvar_par), intent(inout) :: x_guess_gp !initial guess at grid point (in vector form)
    !The new optimal choice will be returned in its place.
    !IMPORTANT: only overwrite it if we actually found something better than the initial guess, otherwise
    !keep it as it is.
    integer, intent(in) :: maxatt !the maximum number of attempts in maximization. Usually this has to be 1, otherwise
    !we spend too much time here. The exception is the first stage of CTS algorithm, where (in the first iteration),
    !we can take quite a lot of attempts (perhaps up to 50) - and this is governed by parameter CTS_first_maxatt in
    !parameter file.

    integer, intent(in) :: gp_ind !this is for debugging only (so that we can write out some
    !additional information for a particular grid point).

    real(dp), intent(in) :: runtime_limit_unix
    logical, intent(in) :: PSO_this_iter

    real(dp), intent(out) :: V_new_gp !New optimum value at the grid point.
    logical, intent(in) :: skip_max !If true, then the maximization step will be skipped and the
    !function returns the same initial guess as was passed to it, and the value obtained
    !using this initial guess.

    real(dp), intent(out) :: x_AARD_gp,x_MARD_gp !The average and maximum absolute relative deviation in policy on a grid point.
    !(+1 in denominator in order to not penalize it excessively when we are close to zero).
    real(dp), dimension(ncvar_par) :: x_ARD_gp !this will not be returned, it is just useful so we do not compute the ARD twice

    real(dp), intent(inout) :: x_AAD_gp !average absolute deviation (change) in policy function.
    !this is used only to adjust rhoend if rhoendadapt == .true. (set in parameter file).

    !Pointer to value function. This will be used to make the old value function V_old
    !accessible to the subroutine used in NAG maximization.
    real(dp), dimension(:,:,:,:,:), pointer :: V_old_pntr

    !Values of state variables as saved in the grid (hence the suffix st).
    real(dp), dimension(I_par) :: a_st
    real(dp), dimension(I_par-1) :: rho_st
    real(dp) :: aex_st
    integer :: g_st_ind

    !Variables used in NAG subroutine E04JYF
    real(dp), dimension(ncvar_par) :: xc !choice variables (M_par x I_par matrix in row form)
    real(dp) :: fc !function value in maximization
    !The following two arrays are arrays of reals and integers which can be used to supply information
    !to the maximized function as an alternative to COMMON block.
    real(dp), dimension(1) :: ruser !Contains no information, only defined for purposes of using a NAG subroutine
    integer, dimension(3) :: iuser !will be used to pass what_binds to some subroutines.
    integer :: ifail, nf
    !Bounds:
    real(dp), dimension(ncvar_par) :: bl
    real(dp), dimension(ncvar_par) :: bu
    real(dp) :: fmax
    !some other things needed by  e04jyf
    !Try changing the dimension of these things and see what happens (there is only a lower bound)
    integer :: LIW = ncvar_par + 2,LW=max(ncvar_par*(ncvar_par-1)/2+12*ncvar_par,13)
    integer, dimension(ncvar_par + 2) :: IW
    real(dp), dimension(max(ncvar_par*(ncvar_par-1)/2+12*ncvar_par,13)) :: W

    real(dp), dimension(ncvar_par) :: bcsf_x !best choice so far
    real(dp) :: bcsf_fmax !best choice so far return (the one given by CPC_ret, so minus one of actual return).

    real(dp) :: IG_value !Value of initial guess (current + continuation)

    !loss - used for reporting statistics on loss resulting from violating constraints
    real(dp) :: loss !Loss for a particular choice passed by CPC_ret common block.
    real(dp), intent(out) :: loss_gp !loss associated with the optimal choice at the grid point.

    integer, dimension(M_par,4) :: j_guess !The initial guess of position of next-period states
    !in the grid, for evry shock realization.

    integer :: att_ind !index for cycling over number of attempts.

    real(dp), dimension(ncvar_par) :: rnd_num

    !_________The following variables are used only in PSO____________
    integer :: ndim, npar
    integer, parameter :: LIOPTS = 100,lopts = 100
    integer :: iopts(LIOPTS)
    real(dp) :: opts(lopts)
    integer :: itt(6)
    real(dp), dimension(ncvar_par) :: bl_pso
    real(dp), dimension(ncvar_par) :: bu_pso
    integer :: inform
    !__________END of PSO variables____________


    !Common block to be accessible by subroutine CPC_ret
    common /cpc_ret_bloc/ a_st,rho_st,aex_st,g_st_ind,j_guess
    !(EFE) : The bounds could actually be prepared even before this subroutine is called,
    !they are the same at every gridpoint.
    !common /cpc_ret_bloc2/ g_lvl_mat, P_onerow,l_max_mat,b_min_mat,b_max_mat,a_min_mat,a_max_mat,rho_min_vec,rho_max_vec
    common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
    !to some subroutines.
    common /commloss/ loss

    !Initialize iuser, ruser (they are not used in any way, just to avoid errors)
    iuser = 0
    ruser = 0.0_dp

    !Much of what follows is similar to what was done elsewhere (get_guess, and LFFC).

    !Get current state from grid:
    a_st(1) = grids%a1_gr(ind_unfold(1))
    a_st(2) = grids%a2_gr(ind_unfold(2))
    rho_st = grids%rho_gr(ind_unfold(3))
    aex_st = grids%aex_gr(ind_unfold(4))
    g_st_ind = ind_unfold(5)

    !Initialize the value of initial guess of position of choice in the grid
    !(we could be quite far from the middle of the grid in some dimensions, so the
    !initial guess is not used the first time that we compute the continuation value).
    j_guess = 0

    !First of all, compute the value of the initial guess (x_guess_gp), using the same value as everywhere else.
    !This will be a lower bounds on value. If the algorithm finds a choice which is worse than the initial guess
    !(due to some error for example), then the returned choice will be the initial guess.

    !Compute the value of initial guess.
    call CPC_ret(ncvar_par,x_guess_gp,IG_value,iuser,ruser)
    !Save the value of the initial guess as the new value function. Initial guess is not
    !updated (i.e., policy function the same)

    !Copy the value of initial guess (we need this now even if attmax = 1 because we are using
    !PSO).
    bcsf_fmax = IG_value

    loss_gp = loss !save the loss at the initial guess (unless something better found later,
    !this will be the point chosen and the associated loss will be the loss at this gp)
    !(loss is passed to this function from CPC_ret using a common block).

    !If we are skipping maximization (Howard's algorithm), then save the value of initial guess as the
    !new value, and return.
    if(skip_max) then
        V_new_gp = -IG_value
        !(CPC_RET returns the negative of value)
        !Change in policy is zero in this case
        x_AARD_gp = 0.0_dp
        x_MARD_gp = 0.0_dp
        return
    end if


    !If we are doing PSO at this iteration, run it, take the optimum found (and check its value). If it is
    !strictly better than the best choice so far (which is just the initial guess), then use it as the new center
    !for permutations. This is pretty much the same as what we do in subroutine sim_series, but it's simpler here
    !because we only save 1 point around which permutations are made, and stage '2a' of the randomization is not
    !done here (it proved pretty useless in simulations). An issue with PSO is that we need to run initialization
    !subroutines, etc. - it's really quite slow. Perhas I could do it only once per VFI iteration per image using
    !the save attribute somehow - but I'm not sure how to do it exactly and if it's possible - investigate this
    !later if PSO seems to be slow... (still the initialization probably isn't where the algorithm spends the most time).
    if(PSO_this_iter) then
        !Initialize this again just in case (to avoid being stuck in an interpolation subroutine).
        j_guess = 0

        ifail = 0
        Call e05zkf('Initialize = E05SAF',iopts,liopts,opts,lopts,ifail)

        !Various parameters may be changed using zkf. See the sample program on website of FL 24.
        !The default values of the following two parameters are E-04. See what different it makes if
        !I set them to E-5 or E-6 in terms of runtime.
        ifail = 0
        Call e05zkf('Local Interior Tolerance = 0.00001',iopts,liopts,opts,lopts,ifail)
        ifail = 0
        Call e05zkf('Local Exterior Tolerance = 0.00001',iopts,liopts,opts,lopts,ifail)

        ndim = 5
        npar = pars%VFI_PSO_par !suggested value is 10*ndim (but perhaps I can try more/less).

        ifail = 1 !silent exit (in fact the documentation suggests this value because ifail = 1
        !every time that it is not guaranteed that global minimum was found - which is most of the times).

        !set bounds
        bl_pso = (1-pars%VFI_PSO_perm)*x_guess_gp
        bu_pso = (1+pars%VFI_PSO_perm)*x_guess_gp

        !If any of the bounds are close to zero, do not even attempt PSO.
        !Here I assume that x_guess is greater than zero (so bl < bu). If not we have other more serious
        !problems, and the algorithm will crash elsewhere anyway.
        if(minval(reshape(abs(bl_pso),[ncvar_par])) < 0.001_dp) then
            goto 666
        end if


        !Call the PSO subroutine
        Call e05saf(ndim,npar,xc,fc,bl_pso,bu_pso,CPC_ret_PSO,E05SXM,iopts,opts, &
            iuser,ruser,itt,inform,ifail)


        !Just in case that the return fc does not corresponds to xc (this can sometimes happen if a nag
        !subroutine encounters an error) - recompute the return.
        call CPC_ret(5,xc,fc,iuser,ruser)


        !Now check if it's better than x_guess_gp. If yes, x_guess_gp. We will then proceed with the local improvement
        !starting at this point in the same way as if PSO was not used
        if(fc < ig_value) then
            x_guess_gp = xc !We don't save this into 'bcsf_x or something like that but we overwrite x_guess_gp).
            bcsf_fmax = fc !best choice so far updated.
        end if

        j_guess = 0 !again, just in case it somehow got overwritten during PSO (It shouldn't happen but if it does
        !it could lead to catastrophic crashes later on)
    666 end if


    !If maxatt == 1 (this will be the case almost always) then we treat this separately (no
    !permutations, etc.) for maximum efficiency
    if(maxatt == 1) then

    ifail = 1 !silent exit

    select case(pars%opt_subroutine)
        case(1) !Quadratic approximation method
            bl = 0.0001_dp !Set lower bounds
            xc = x_guess_gp !get initial guess

            !Notice the adaptation of rhoend. 0.1 multiple of the average change in policy
            !function (in absolute value) across the choice variables, in the previous iteration of VFI.
            !We might change it to 0.01 multiple in the future to see what it changes.
            call e04jcf(CPC_ret2,ncvar_par,ncvar_par*2+1,xc,bl,pars%x_bu,0.05_dp,&
            min(pars%rhoend,0.1_dp*x_AAD_gp),e04jcp,pars%VFI_opt_maxiter,fmax,nf,iuser,ruser,ifail)
        case(2) !Newton-like method - this tends to lead to crashes on HPC for some reason.
            bl = 0.0001_dp
            bu = pars%x_bu
            xc = x_guess_gp
            call e04jyf(ncvar_par,0,CPC_ret,bl,bu,xc,fmax,iw,liw,w,lw,iuser,ruser,ifail)
        case default
            write(*,*) 'Error in subroutine VFI_solve_gp - wrong value of optim_subroutine'
            error stop
    end select
    !We assume that fmax is the value of CPC_ret associated with choice xc. This should be
    !the case always according to documentation of the NAG subroutines e04jcf,e04jyf, but it might be safer
    !to recompute the value of xc just in case something bad happened.

    !Check whether the optimum found by one of the subroutines is better than the
    !value of the initial guess. If yes, update the initial guess and use its
    !value as the new value. Otherwise do not update the initial guess and
    !use its value.

    if(fmax<bcsf_fmax) then !(optimum found is better than init. guess (or what PSO found), remember we work with (-1) multiples)
        !BEFORE overwriting the initial guess x_guess_gp with the new one, compute the difference.
        !This will be used one level higher to compute various statistics for convergence/debugging.
        x_ARD_gp = abs((xc - x_guess_gp)/(abs(x_guess_gp) + 1.0_dp))

        x_AARD_gp = sum(x_ard_gp)/real(ncvar_par,dp)
        x_MARD_gp = maxval(x_ARD_gp)

    !IF we are adjusting rho, also compute the change in absolute value of policy (levels), averaged
    !across the choice variables. We need to work with levels here. This could sometimes be zero (in case
    !we found the exact same point as before) - in that case make sure that rhoend is not less than, say
    !E-18, which corresponds to x_AAT_gp = E-17. We never want rhoend = 0.
    x_AAD_gp = max(sum(abs(xc - x_guess_gp))/real(ncvar_par,dp),1.0E-17_dp);

        !(note): Here we only perform one maximization, so we just check whether the function value is better
        !than the value of the initial guess. In the future, we could try several permutations of the initial guess,
        !and/or more optimization subroutines.
        x_guess_gp = xc
        V_new_gp = -fmax

        !WE also need to compute the loss associated with this point. To do this use CPC_ret, which
        !will save loss. This is for debugging purposes (and for detecting divergence issues).
        call CPC_ret(ncvar_par,x_guess_gp,fmax,iuser,ruser)
        loss_gp = loss !(loss is obtained using common block)
    else !Use initial guess instead (and don't update initial guess)
        !We already have the value of initial guess in IG_value, and loss_gp was also saved there (it was
        !not updated since), and the initial guess x_guess_gp was not updated (so it will be returned).
        !We just need to update the value.
        V_new_gp = -IG_value
        !Change in policy is zero in this case
        x_AARD_gp = 0.0_dp
        x_MARD_gp = 0.0_dp

        !In this case x_ARD_gp is not updated (we use the same as the last iteration).
    end if

    else !Here is the branch of the program for if maxatt > 1

    !Copy the initial guess as the best choice so far, and copy the value. Then
    !maximize attmax number of times.
    bcsf_x = x_guess_gp !x_guess_gp should contain the initial guess at his grid point and should not
    !have been overwritten at this stage.

    !This is now redudnant, we save this every time due to introduction of PSO
    !bcsf_fmax = IG_value

    do att_ind = 1,maxatt
    !Here permutate the initial guess unless this is the first attempt, in which case - don't.
    if(att_ind > 1) then
    !rnd_num contains 5 numbers to be used for permutation. of xc_init.
        call random_number(rnd_num)
        rnd_num = 1.0_dp + (rnd_num - 0.5_dp)*pars%VFI_perm !VFI_perm gives percentage deviation bounds. For example
        !if this is 0.1, the deviation is plus minus 10 percent.

        xc = x_guess_gp * rnd_num

        !Also, when we are permutating, it makes little sense to use the adjusted rho based
        !(which is useful if we are in the correct local optimum - but if we are exploring
        !regions further from the current local optimum we don't need to be as accurate - a rhoend
        !of around 1E-06 which is the default should suffice. Perhaps even E-05 for extra speed.
        !I temporarily set x_AAD_gp = 1, which will force the optimization algorithm to use the
        !default value of rhoend. x_AAD_gp will be recomputed at the end. We could even
        !handle this case separately below and call the subroutine with a different rhoend
        !(E-04-E-05 for extra speed - I might try it later).
        x_AAD_gp = 1.0_dp
    else
        xc = x_guess_gp
    end if

    ifail = 1 !silent exit
    select case(pars%opt_subroutine)
        case(1) !Quadratic approximation method
            bl = 0.0001_dp !Set lower bounds

            call e04jcf(CPC_ret2,ncvar_par,ncvar_par*2+1,xc,bl,pars%x_bu,0.05_dp,&
            min(pars%rhoend,0.1_dp*x_AAD_gp),e04jcp,pars%VFI_opt_maxiter,fmax,nf,iuser,ruser,ifail)
        case(2) !Newton-like method - this tends to lead to crashes on HPC for some reason.
            bl = 0.0001_dp
            bu = pars%x_bu
            call e04jyf(ncvar_par,0,CPC_ret,bl,bu,xc,fmax,iw,liw,w,lw,iuser,ruser,ifail)
        case default
            write(*,*) 'Error in subroutine VFI_solve_gp - wrong value of optim_subroutine'
            error stop
    end select
    !We assume that fmax is the value of CPC_ret associated with choice xc. This should be
    !the case always according to documentation of the NAG subroutines e04jcf,e04jyf, but it might be safer
    !to recompute the value of xc just in case something bad happened.

    !Check whether the optimum found by one of the subroutines is better than the best choice found
    !so far. IF yes, overwrite the previously found best choice
    if(fmax < bcsf_fmax) then
        bcsf_fmax = fmax
        bcsf_x = xc
    end if

    !If the runtime limit is reached, quit now (no more attempts).
    if(time_real_acc() > runtime_limit_unix) then
        exit
    end if

    end do

    if(bcsf_fmax<IG_value) then !(optimum found is better than init. guess, remember we work with (-1) multiples)
        !Update initial guess (have to reshape it into incomplete cons. matrix first)

        !BEFORE overwriting the initial guess x_guess_gp with the new one, compute the difference.
        !This will be used one level higher to compute various statistics for convergence/debugging.
        x_ARD_gp = abs((bcsf_x - x_guess_gp)/(abs(x_guess_gp) + 1.0_dp))

        x_AARD_gp = sum(x_ard_gp)/real(ncvar_par,dp)
        x_MARD_gp = maxval(x_ARD_gp)


    !IF we are adjusting rho, also compute the change in absolute value of policy (levels), averaged
    !across the choice variables. We need to work with levels here. This could sometimes be zero (in case
    !we found the exact same point as before) - in that case make sure that rhoend is not less than, say
    !E-18, which corresponds to x_AAT_gp = E-17. We never want rhoend = 0.
    x_AAD_gp = max(sum(abs(bcsf_x - x_guess_gp))/real(ncvar_par,dp),1.0E-17_dp);

        !(note): Here we only perform one maximization, so we just check whether the function value is better
        !than the value of the initial guess. In the future, we could try several permutations of the initial guess,
        !and/or more optimization subroutines.
        x_guess_gp = bcsf_x
        V_new_gp = -bcsf_fmax

        !WE also need to compute the loss associated with this point. To do this use CPC_ret, which
        !will save loss. This is for debugging purposes (and for detecting divergence issues).
        call CPC_ret(ncvar_par,x_guess_gp,fmax,iuser,ruser)
        loss_gp = loss !(loss is obtained using common block)
    else !Use initial guess instead (and don't update initial guess)
        !We already have the value of initial guess in IG_value, and loss_gp was also saved there (it was
        !not updated since), and the initial guess x_guess_gp was not updated (so it will be returned).
        !We just need to update the value.
        V_new_gp = -IG_value
        !Change in policy is zero in this case
        x_AARD_gp = 0.0_dp
        x_MARD_gp = 0.0_dp

        !If we did not find a better choice than previously, then do not upate the X_AAD_gp.
    end if


    end if !End of the branch of program with more than one initial guess attempt (maxatt > 1)



end subroutine VFI_solve_gp


!Subroutine solve_FP_problem solves the first period problem (t=0), given a value function
!for periods t>=1 and initial conditions.
!
!It is better to have one subroutine for solving the FP problem and one for solving
!the subsequent periods because the dimensionality of choice variables is different
!and there are a lot of other small changes which would make it confusing if we
!put all simulation into one subroutine).
subroutine solve_FP_problem(V_old,b_init,b_ex_init,g_ind_init,pars,SS,FP_result)
    Use nag_library, Only: e04jyf, nag_wp, e04jcf, e04jcp,e05jaf,e05jbf,E05JBK
    real(dp), dimension(:,:,:,:,:), intent(in) :: V_old
    real(dp), dimension(1,I_par), intent(in) :: b_init !initial debt (it can be permutated, hence don't use pars%b_init)
    real(dp), intent(in) :: b_ex_init !initial external debt
    integer, intent(in) :: g_ind_init !initial shock realization
    type(par), intent(in) :: pars
    type(ss_type), intent(in) :: SS !autarky steady state (for getting an initial guess).
    type(FP_result_type), intent(out) :: FP_result

    !Variables used in NAG subroutine E04JYF
    real(dp), dimension(3) :: xc,xc2,xc_tmp, xc_init !choice variables

    real(dp) :: fc !function value in maximization
    !The following two arrays are arrays of reals and integers which can be used to supply information
    !to the maximized function as an alternative to COMMON block.
    real(dp), dimension(1) :: ruser !Contains no information, only defined for purposes of using a NAG subroutine
    integer, dimension(3) :: iuser !will be used to pass what_binds to some subroutines.
    integer :: ifail
    !Bounds:
    real(dp), dimension(3) :: bl = 0.0001_dp,bu

    integer :: att_ind
    real(dp),dimension(3) :: rnd_num

    !TMP variables needed for testing E04jcf (remove/comment if I don't use that subroutine
    !in the end)
    integer :: nf
    !Local
    logical :: c_im_found
    real(dp) :: IG_value !Value of initial guess (current + continuation)

    !Best choice so far
    real(dp), dimension(3) :: bcsf_x !choice
    real(dp) :: bcsf_fc !value (as returned by subroutine, so we are looking for the smallest possible one).

    !Copies of data which we will pass to CPC subroutine using common block (need
    !a copy because these object are passed as arguments to the subroutine
    !so they can't be in a common block).
    real(dp), dimension(1,I_par) :: b_init_cp
    real(dp) :: b_ex_init_cp
    integer :: g_ind_init_cp

    integer, dimension(4) :: j_guess !The guess of initial position of a point on the grid. This will be common
    !to subroutines solve_FP_problem, and the subroutine CPC_ret_FP. It will be intialized at zero in subroutine
    !solve_FP_problem, then the new initial guess will be updated every time that CPC_ret is called by
    !a maximization subroutine (then it calls subroutine interp_V, which returns the interpolated value, and
    !updates the guess).

    real(dp) :: tmp_debug !can delete this later

    real(dp), dimension(1,1) :: rho_prime,aex_prime,b_ex_prime !temporary variables (because dimension in
    !FP_result and in get_av_fp is not consistent

    !Try changing the dimension of these things and see what happens (there is only a lower bound)
    integer :: LIW = 5,LW=max(3*(3-1)/2+12*3,13)
    integer, dimension(5) :: IW
    real(dp), dimension(max(3*(3-1)/2+12*3,13)) :: W

    !The states are passed to CPC_ret_fp function using a commong block
    !(they may be permutated)
    common /cpc_ret_fp_bloc/ b_init_cp, b_ex_init_cp ,g_ind_init_cp,j_guess

    !Create copies of inputs which need to be made available to subroutine
    !called by the maximization subroutine using a common block.
    b_init_cp = b_init
    b_ex_init_cp = b_ex_init
    g_ind_init_cp = g_ind_init

    !initialize iuser, ruser, nf
    iuser = 0
    ruser = 0.0_dp
    nf = 0

    !Save the actual index in FP_result (there can be randomization over this)
    FP_result%s_ind_init = g_ind_init_cp

    !Generate bounds ul and ub
    bl = 0.00001_dp !Lower bounds

    !With two countries, there are now three free choice variables (consumption of country 1,
    !consumption of country 2, and labour supply of country 1).

    !We start with an initial guess which corresponds to autarky steady-state and then we
    !permutate around this randomly. Save this into x_init.
    xc_init(1) = SS%l(1,1) !First element is labour supply of country 1
    xc_init(2) = SS%c(1,1) !Second element is consumption of country 1
    xc_init(3) = SS%c(1,2) !Second element is consumption of country 2




    !Initialize best choice found so far
    bcsf_x = xc_init

    !IF we are only interested in getting a SS with trasnfers (then the autarky SS guess) is the
    !the correct choice in this regard. Skip the rest of the subroutine.
    if(pars%gen_transfers_ss) then
        goto 545
    end if

    !And the value of initial guess (bcsf_fc contains the negative of return, so lower values are better!)
    call CPC_ret_FP(3,xc_init,bcsf_fc,iuser,ruser,nf)

    !debug - save the value of initial guess for later comparison with actually found optimum
    tmp_debug = bcsf_fc

    !Get bounds for maximization. Will use the same bounds as for the optimization in periods t>1,
    !which are already saved in parameter file. Just need to recover the appropriate elements.
    if(I_par /= 2) then
        write(*,*) 'Error: subroutine solve_FP_problem needs to be generalized for I_par /= 2'
        error stop
    end if
    bu(1) = pars%x_bu(1) !bound for labour supply in country 1
    bu(2) = pars%x_bu(1+M_par) !bound for consumption in country 1
    if(M_par > 1) then !bound for consumption in country 2
        bu(3) = pars%x_bu(1+2*M_par)
    else
        bu(3) = bu(2) !If deterministic env, then cons in country 2 isn't part of xc in period t>1.
        !just use the consumption in country 1 bounds instead
    end if

    !Attempt (possibly many times) maximizations. If pars%perm_eps = 0, then we are solving the initial guess
    !problem only once, and we can afford more attempts.
    do att_ind = 1,pars%sim_attmax1
        !If this is the first attempt, do not permutate. Otherwise permutate around xc_init
        if(att_ind > 1) then
            !rnd_num contains 3 numbers to be used for permutation. of xc_init.
            call random_number(rnd_num)
            !We want to transform this into interval (1-pars%sim_perm1,1+pars%sim_perm1)
            rnd_num = rnd_num * pars%sim_perm1 * 2.0_dp + 1.0_dp - pars%sim_perm1

            xc = xc_init * rnd_num
        else
            xc = xc_init
        end if

        !Initialize the guess of position in the grid. Two alternatives are attractive. Either we can start in the
        !middle of the grids (that is where next-period states should be close to unless permutations are big,
        !or we can just take an initial guess equal to zero (so that the first call of CPC_ret_FP, the guess
        !is not used, which is more efficient if the guess is bad). It doesn't matter very much.
        j_guess(1) = pars%N_a / 2
        j_guess(2) = pars%N_a / 2
        j_guess(3) = pars%N_rho / 2
        j_guess(4) = pars%N_aex / 2

        !Perform maximization
        ifail = 1 !silent exit
        !ifail = -1 !noisy exit
        select case(pars%opt_subroutine)
            case(1) !Quadratic approximation method
                bl = 0.0001_dp !Set lower bounds
!                call e04jcf(CPC_ret_FP,3,7,xc,bl,bu,0.05_dp,&
!                    pars%rhoend*0.01_dp,e04jcp,pars%sim_maxiter,fc,nf,iuser,ruser,ifail)

                call e04jcf(CPC_ret_FP,3,7,xc,bl,bu,0.05_dp,&
                    pars%sim_rhoend,e04jcp,pars%sim_maxiter,fc,nf,iuser,ruser,ifail)
            case(2) !Newton-like method - this tends to lead to crashes on HPC for some reason.
                bl = 0.0001_dp
            call e04jyf(3,0,CPC_ret_FP2,bl,bu,xc,fc,iw,liw,w,lw,iuser,ruser,ifail)
                !Check subroutine vfi_solve_gp for details on generalization, it's used there.
            case default
                write(*,*) 'Error in subroutine Solve_FP_problem - wrong value of optim_subroutine'
                error stop
        end select

        !If the function value is better (lower because minimization) than the best one so far,
        !save the new best value, and the best guess
        if(fc < bcsf_fc) then
            bcsf_fc = fc
            bcsf_x = xc
        end if
    end do

    !debug - check out how the best choice varies with different values of parameters (and perhaps
    !optimization subroutines)


    !now bscf_x contains the best choice found (in vector form). We need to recover
    !all the remaining variables using subroutine get_av_fp, and save these in appropriate
    !elements of FP_result. Then we are done.

545 FP_result%s_ind_init = g_ind_init
    !From the optimal choice, get all the remaining variables and save them in appropriate elements of FP_result
    call get_av_fp(pars,bcsf_x,b_init,b_ex_init,g_ind_init,ifail,&
        FP_result%c,FP_result%l,FP_result%b_prime,b_ex_prime,FP_result%a_prime,rho_prime,aex_prime,.false.)

    !Get copies of the inputs with different dimension
    FP_result%b_ex_prime = b_ex_prime(1,1)
    FP_result%rho_prime = rho_prime(1,1)
    FP_result%aex_prime = aex_prime(1,1)

    !Also get marginal utilities of consumption and labour in both countries. Could be useful for debugging
    call util_c(FP_result%c,FP_result%u_c,pars%A_par,pars%sigma_par)
    call util_l(FP_result%l,FP_result%u_l,pars%B_par,pars%Gamma_par)

end subroutine solve_FP_problem

!Subroutine CPC_ret_FP is like CPC_ret (it gives the current and continuation return of a choice).
!The difference is that this is for first period, in which the constraints are slightly different,
!and the choice variables are different.
subroutine CPC_ret_FP(n,xc,fc,iuser,ruser,inform)
        real (dp), intent (out) :: fc
        integer, intent (in) :: n
        real (dp), intent (inout) :: ruser(*)
        real (dp), intent (in) :: xc(n)
        integer, intent (inout) :: iuser(*) !In this instance, it is used to pass what_binds
        integer, intent(out) :: inform

        !initial states
        real(dp), dimension(1,I_par) :: b_init_cp
        real(dp) :: b_ex_init_cp
        integer :: g_ind_init_cp

        real(dp), dimension(:,:,:,:,:), pointer :: V_old_pntr !pointer to V_old, passed using common block
        real(dp), dimension(:), pointer :: a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr !pointers so all grids

        !pointer pars_pntr points to parameters type (it is more convenient to pass it using a common
        !block than create copies of everything we need and put them in a long common block). We
        !just need to be careful to not make any assignment statements because we would then modify
        !parameters elsewhere in the program!!!
        type(par), pointer :: pars_pntr

        integer :: fail
        real(dp), dimension(1,I_par) :: c,l,a_pr
        real(dp), dimension(1,1) :: rho_pr, aex_pr !dimension for minimal changes in code relative to CPC_ret
        !fail,c,l,a_pr,rho_pr,aex_pr,get_cl_only

        real(dp) :: loss
        logical :: constr_ok
        real(dp) :: conval !continuation value

        !asset holdings without MU adjustment.
        real(dp), dimension(1,I_par) :: b_pr
        real(dp), dimension(1,1) :: b_ex_pr

        !Utility in both countries
        real(dp), dimension(1,I_par) :: U

        integer, dimension(4) :: j_guess !guess of position of the point on the grid.

        !Current states, value function, parameters, and grids are passed to this function
        !using common block (and sometimes pointers too)
        common /cpc_ret_fp_bloc/ b_init_cp, b_ex_init_cp ,g_ind_init_cp,j_guess
        common /V_old_pointer/ V_old_pntr
        common /pars_pointer/ pars_pntr
        common /grids/ a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr

        !The guess j_guess should be in the common block.

        real(dp), dimension(M_par,1), parameter ::  ones = 1.0_dp

        !Negative value of inform leads to termination of maximization procedure. Never use this, it leads to crashes.
        inform = 0
        loss = 0.0_dp

        !Test subroutine get_av_fp(pars,x,b_init,b_ex_init,g_st_ind,fail,c,l,a_pr,rho_pr,aex_pr,get_cl_only)
        call get_av_fp(pars_pntr,xc,b_init_cp,b_ex_init_cp,g_ind_init_cp,fail,c,l,b_pr,b_ex_pr,a_pr,rho_pr,aex_pr,.false.)

        if(fail > 0) then
            fc = max_loss !maximum loss. This introduces a significant discontinuity, need to be careful about implications elsewhere
            return
        end if

        !If we got here, it means that the point is reasonable at least in the sense of satisfying the non-negativity constraint.
        !Check all remaining constraints and compute the loss (unadjusted).

        !We still want to use subroutine check_const, which works with state-contingent plans). So we stack
        !each row vector on top of itself M_par times to get a state-contingent plan, and then divide the
        !resulting loss by M_par to get the correct loss value.

        call check_constr(pars_pntr,matmul(ones,c),matmul(ones,l),matmul(ones,a_pr)&
        ,matmul(ones,rho_pr),matmul(ones,aex_pr),constr_ok,loss,.false.,.false.)
        loss = loss/real(M_par,dp) !normalize so the loss is penalized the same way as in periods 1,...
        !(it doesn't really matter, in the first period we should never choose a point for which loss is positive)

        !Now that we know consumption and labour plan, we can compute the current return. There is no expectation
        !unlike in period t>=1 problem, so do this manually.
        !Get utility of both countries
        call util(c,l,u,pars_pntr%A_par,pars_pntr%B_par,pars_pntr%sigma_par,pars_pntr%Gamma_par)

        !Current return
        fc = sum(pars_pntr%mass*pars_pntr%alpha*U)

        !Now we have current return in fc, and loss in loss. We just need to compute the continuation return.

        !If the choice is not strictly feasible, force it to be within bounds (extrapolation
        !is very bad, usually leads to divergence). So effectively if the algorithm chooses a point
        !outside of the grid boundaries, we use neares neighbour extapolation (sort of) and add
        !a loss. This is much more stable (otherwise in some corners of the grid the loss if we
        !use extrapolation is compounded and the program diverges). This puts a lower bound on it and it is
        !fine as long as the region where this happens was not optimally reached in the first place (need
        !to experiment with loss).
        if(.not. constr_ok) then
            a_pr=max(pars_pntr%a_min,a_pr)
            a_pr=min(pars_pntr%a_max,a_pr)
            rho_pr=max(pars_pntr%rho_min,rho_pr)
            rho_pr=min(pars_pntr%rho_max,rho_pr)
            aex_pr=max(pars_pntr%aex_min,aex_pr)
            aex_pr=min(pars_pntr%aex_max,aex_pr)
        end if

        !For testing purposes, it is sometimes useful to try quadrilinear interpolation in the first period
        if(pars_pntr%sim_FP_quadrilin) then
               call interp_V(V_old_pntr,[a_pr(1,1),a_pr(1,2),rho_pr(1,1),aex_pr(1,1)],j_guess&
                ,min(g_ind_init_cp,pars_pntr%M_grid),conval,a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr,size(a1_gr_pntr),&
                size(rho_gr_pntr),size(aex_gr_pntr))
        else

        !Which interpolation function we use depends on the value of VFI_interp_mode (set in parameter file).
        select case(pars_pntr%VFI_interpolation_mode)
            case(1) !Quadrilinear interpolation
               call interp_V(V_old_pntr,[a_pr(1,1),a_pr(1,2),rho_pr(1,1),aex_pr(1,1)],j_guess&
                ,min(g_ind_init_cp,pars_pntr%M_grid),conval,a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr,size(a1_gr_pntr),&
                size(rho_gr_pntr),size(aex_gr_pntr))
            case(2) !Local Shepard's interpolation
               call interp_V_shep(V_old_pntr,[a_pr(1,1),a_pr(1,2),rho_pr(1,1),aex_pr(1,1)],j_guess&
                ,min(g_ind_init_cp,pars_pntr%M_grid),conval,a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr,size(a1_gr_pntr),&
                size(rho_gr_pntr),size(aex_gr_pntr),pars_pntr%shep_n,pars_pntr%shep_norm,pars_pntr%shep_norm_val,pars_pntr%shep_floor)
            case(3)
               call interp_V_spline(V_old_pntr,[a_pr(1,1),a_pr(1,2),rho_pr(1,1),aex_pr(1,1)],&
                min(g_ind_init_cp,pars_pntr%M_grid),conval,a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr,size(a1_gr_pntr),&
                size(rho_gr_pntr),size(aex_gr_pntr),pars_pntr%bspline_k_a1,pars_pntr%bspline_k_a2,pars_pntr%bspline_k_rho,pars_pntr%bspline_k_aex)
        end select

        endif

        !Now add continuation value (discounted) to fc, penalize by subtracting loss,
        !and multiply by (-1) because we are working with a minimization subroutine, so we want to 'minimize -f'
        fc = (fc + pars_pntr%beta*conval - loss) * (-1.0_dp)

end subroutine CPC_ret_FP


!Subroutine CPC_ret_FP is a copy of CPC_ret_FP, but with a different interface (so that a different
!NAG subroutine can use it). This could be done nices (less code) but it is low priority right now.
subroutine CPC_ret_FP2(n,xc,fc,iuser,ruser)
        real (dp), intent (out) :: fc
        integer, intent (in) :: n
        real (dp), intent (inout) :: ruser(*)
        real (dp), intent (in) :: xc(n)
        integer, intent (inout) :: iuser(*) !In this instance, it is used to pass what_binds

        !initial states
        real(dp), dimension(1,I_par) :: b_init_cp
        real(dp) :: b_ex_init_cp
        integer :: g_ind_init_cp

        real(dp), dimension(:,:,:,:,:), pointer :: V_old_pntr !pointer to V_old, passed using common block
        real(dp), dimension(:), pointer :: a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr !pointers so all grids

        !pointer pars_pntr points to parameters type (it is more convenient to pass it using a common
        !block than create copies of everything we need and put them in a long common block). We
        !just need to be careful to not make any assignment statements because we would then modify
        !parameters elsewhere in the program!!!
        type(par), pointer :: pars_pntr

        integer :: fail
        real(dp), dimension(1,I_par) :: c,l,a_pr
        real(dp), dimension(1,1) :: rho_pr, aex_pr !dimension for minimal changes in code relative to CPC_ret
        !fail,c,l,a_pr,rho_pr,aex_pr,get_cl_only

        real(dp) :: loss
        logical :: constr_ok
        real(dp) :: conval !continuation value

        !asset holdings without MU adjustment.
        real(dp), dimension(1,I_par) :: b_pr
        real(dp), dimension(1,1) :: b_ex_pr

        !Utility in both countries
        real(dp), dimension(1,I_par) :: U

        integer, dimension(4) :: j_guess !guess of position of the point on the grid.

        !Current states, value function, parameters, and grids are passed to this function
        !using common block (and sometimes pointers too)
        common /cpc_ret_fp_bloc/ b_init_cp, b_ex_init_cp ,g_ind_init_cp,j_guess
        common /V_old_pointer/ V_old_pntr
        common /pars_pointer/ pars_pntr
        common /grids/ a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr

        !The guess j_guess should be in the common block.

        real(dp), dimension(M_par,1), parameter ::  ones = 1.0_dp

        loss = 0.0_dp

        !Test subroutine get_av_fp(pars,x,b_init,b_ex_init,g_st_ind,fail,c,l,a_pr,rho_pr,aex_pr,get_cl_only)
        call get_av_fp(pars_pntr,xc,b_init_cp,b_ex_init_cp,g_ind_init_cp,fail,c,l,b_pr,b_ex_pr,a_pr,rho_pr,aex_pr,.false.)

        if(fail > 0) then
            fc = max_loss !maximum loss. This introduces a significant discontinuity, need to be careful about implications elsewhere
            return
        end if

        !If we got here, it means that the point is reasonable at least in the sense of satisfying the non-negativity constraint.
        !Check all remaining constraints and compute the loss (unadjusted).

        !We still want to use subroutine check_const, which works with state-contingent plans). So we stack
        !each row vector on top of itself M_par times to get a state-contingent plan, and then divide the
        !resulting loss by M_par to get the correct loss value.

        call check_constr(pars_pntr,matmul(ones,c),matmul(ones,l),matmul(ones,a_pr)&
        ,matmul(ones,rho_pr),matmul(ones,aex_pr),constr_ok,loss,.false.,.false.)
        loss = loss/real(M_par,dp) !normalize so the loss is penalized the same way as in periods 1,...
        !(it doesn't really matter, in the first period we should never choose a point for which loss is positive)

        !Now that we know consumption and labour plan, we can compute the current return. There is no expectation
        !unlike in period t>=1 problem, so do this manually.
        !Get utility of both countries
        call util(c,l,u,pars_pntr%A_par,pars_pntr%B_par,pars_pntr%sigma_par,pars_pntr%Gamma_par)

        !Current return
        fc = sum(pars_pntr%mass*pars_pntr%alpha*U)

        !Now we have current return in fc, and loss in loss. We just need to compute the continuation return.

        !If the choice is not strictly feasible, force it to be within bounds (extrapolation
        !is very bad, usually leads to divergence). So effectively if the algorithm chooses a point
        !outside of the grid boundaries, we use neares neighbour extapolation (sort of) and add
        !a loss. This is much more stable (otherwise in some corners of the grid the loss if we
        !use extrapolation is compounded and the program diverges). This puts a lower bound on it and it is
        !fine as long as the region where this happens was not optimally reached in the first place (need
        !to experiment with loss).
        if(.not. constr_ok) then
            a_pr=max(pars_pntr%a_min,a_pr)
            a_pr=min(pars_pntr%a_max,a_pr)
            rho_pr=max(pars_pntr%rho_min,rho_pr)
            rho_pr=min(pars_pntr%rho_max,rho_pr)
            aex_pr=max(pars_pntr%aex_min,aex_pr)
            aex_pr=min(pars_pntr%aex_max,aex_pr)
        end if

        !Get continuation value using interpolation. This version does not use the initial guess.
!        call interp_V(V_old_pntr,[a_pr(1,1),a_pr(1,2),rho_pr(1,1),aex_pr(1,1)]&
!            ,min(g_ind_init_cp,pars_pntr%M_grid),conval,a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr,size(a1_gr_pntr),&
!            size(rho_gr_pntr),size(aex_gr_pntr))


        !Which interpolation function we use depends on the value of VFI_interp_mode (set in parameter file).
        select case(pars_pntr%VFI_interpolation_mode)
            case(1) !Quadrilinear interpolation
               call interp_V(V_old_pntr,[a_pr(1,1),a_pr(1,2),rho_pr(1,1),aex_pr(1,1)],j_guess&
                ,min(g_ind_init_cp,pars_pntr%M_grid),conval,a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr,size(a1_gr_pntr),&
                size(rho_gr_pntr),size(aex_gr_pntr))
            case(2) !Local Shepard's interpolation
               call interp_V_shep(V_old_pntr,[a_pr(1,1),a_pr(1,2),rho_pr(1,1),aex_pr(1,1)],j_guess&
                ,min(g_ind_init_cp,pars_pntr%M_grid),conval,a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr,size(a1_gr_pntr),&
                size(rho_gr_pntr),size(aex_gr_pntr),pars_pntr%shep_n,pars_pntr%shep_norm,pars_pntr%shep_norm_val,pars_pntr%shep_floor)
            case(3)
               call interp_V_spline(V_old_pntr,[a_pr(1,1),a_pr(1,2),rho_pr(1,1),aex_pr(1,1)],&
                min(g_ind_init_cp,pars_pntr%M_grid),conval,a1_gr_pntr,a2_gr_pntr,rho_gr_pntr,aex_gr_pntr,size(a1_gr_pntr),&
                size(rho_gr_pntr),size(aex_gr_pntr),pars_pntr%bspline_k_a1,pars_pntr%bspline_k_a2,pars_pntr%bspline_k_rho,pars_pntr%bspline_k_aex)
        end select

        !Now add continuation value (discounted) to fc, penalize by subtracting loss,
        !and multiply by (-1) because we are working with a minimization subroutine, so we want to 'minimize -f'
        fc = (fc + pars_pntr%beta*conval - loss) * (-1.0_dp)

end subroutine CPC_ret_FP2

subroutine shocks_gen_and_save(folder_name,pars)
    character*(*), intent(in) :: folder_name
    type(par), intent(in) :: pars

    integer, allocatable :: s_ind_all(:,:)
    real(dp) :: k,randnum
    character*256 :: shell_cmd,file_path

    integer :: n_ind,t_ind

    real(dp), dimension(pars%shocks_T) :: s_ind_avg
    real(dp) :: s_ind_avg2

    !allocate memory (first index is time, second index is series)
    allocate(s_ind_all(pars%shocks_T,pars%shocks_N))

    s_ind_avg = 0.0_dp

    !generate shocks (using the transition matrix P or P_sim, depending on the preference).

    do n_ind = 1,pars%shocks_N
    !In the 'first period', randomize over initial shock (unless initial shock given). Then pick the appropriate element
    !of FP_result.
    if(pars%s_init_index == -1) then
        if(M_par == 1) then
            s_ind_all(1,n_ind) = 1
        elseif(M_par == 2) then
            k = pars%P(M_par,1)/(1.0_dp + pars%P(M_par,1) - pars%P(1,1))
            !k is the asymptotic (stationary) probability of state 1. So draw a number from U(0,1),
            !and if it's less than k, the initial state will be set to 1, otherwise it will be set to two.
            call random_number(randnum)
            if(randnum < k) then
                s_ind_all(1,n_ind) = 1
            else
                s_ind_all(1,n_ind) = 2
            end if
        else
            write(*,*) 'Error: need to generalize subroutine shocks_gen_and_solve for M_par>2.'
            error stop
        end if
    else
        s_ind_all(1,n_ind) = pars%s_init_index
    end if

    s_ind_avg(1) = s_ind_avg(1) + s_ind_all(1,n_ind)

    if(pars%use_P_sim) then
        s_ind_all(2:pars%shocks_T,n_ind) = markov(s_ind_all(1,n_ind),pars%P_sim,pars%shocks_T-1)
    else !this is the default option which should be used all the time outside of debugging
        s_ind_all(2:pars%shocks_T,n_ind) = markov(s_ind_all(1,n_ind),pars%P,pars%shocks_T-1)
    end if

    do t_ind = 2,pars%shocks_T
        s_ind_avg(t_ind) = s_ind_avg(t_ind) + s_ind_all(t_ind,n_ind)
    end do

    end do


    !Average shock realization every period

    s_ind_avg = s_ind_avg / real(pars%shocks_N,dp)

    !Average across 100 periods

    s_ind_avg2 = sum(s_ind_avg(1:min(100,pars%shocks_T)))/ real(min(100,pars%shocks_T),dp)


    write(*,*) 'Average value of shock_realization every period for the first 100 periods'
    do t_ind = 1,min(100,pars%shocks_T)
        write(*,*) 't_ind = ,',t_ind,' s_ind_avg(t_ind) = ',s_ind_avg(t_ind)
    end do


     write(*,*) 'Average across the first 100 periods = ',s_ind_avg2

    !Now save the shocks:
    file_path = 'results/'//trim(folder_name)//'/shocks.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) s_ind_all
    close(unit = 20)

    deallocate(s_ind_all)
end subroutine shocks_gen_and_save

subroutine load_shocks(folder_name,pars,s_ind_loaded)
    character*(*), intent(in) :: folder_name
    type(par), intent(in) :: pars
    integer, intent(inout) :: s_ind_loaded(:,:)

    character*256 :: file_path


    !No error checking here - if the file is not existent or the path is wrong,
    !the program will crash. The same happens if the wrong shocks_N and shocks_t are set in the parameter file
    !(they must be the same as those used to generate the shocks).
    write(*,*) 'Loading shocks from file: results/',trim(pars%shocks_path)


    file_path = 'results/'//trim(pars%shocks_path)
    open(unit=20, file = file_path, status = 'old')
    read(unit=20, fmt=*)  s_ind_loaded
    close(unit = 20)

    write(*,*) 'Shocks were successfully loaded'

end subroutine load_shocks


!Subroutine sim_result_alloc allocates memory to allocatable arrays of
!sim_result_type. It contains no checks of whether the arrays are already allocated
!so it needs to be used carefully, otherwise it will lead to crashes, etc.
subroutine sim_result_alloc(sim_result,pars,wtd)
    type(sim_result_type), intent(inout) :: sim_result
    type(par), intent(in) :: pars
    integer, intent(in) :: wtd !what to do (1 = allocate, 0 = deallocate)

    if(wtd == 1) then
        !Allocate memory to simulated series
        allocate(sim_result%c_sim(pars%T_sim,I_par))
        allocate(sim_result%l_sim(pars%T_sim,I_par))
        allocate(sim_result%b_sim(pars%T_sim,I_par))
        allocate(sim_result%b_prime_sim(pars%T_sim,I_par))
        allocate(sim_result%a_prime_sim(pars%T_sim,I_par))
        allocate(sim_result%u_c_sim(pars%T_sim,I_par))
        allocate(sim_result%u_l_sim(pars%T_sim,I_par))
        allocate(sim_result%U_sim(pars%T_sim,I_par))
        allocate(sim_result%rho_prime_sim(pars%T_sim))
        allocate(sim_result%aex_prime_sim(pars%T_sim))
        allocate(sim_result%b_ex_sim(pars%T_sim))
        allocate(sim_result%b_ex_prime_sim(pars%T_sim))
        allocate(sim_result%R_sim(pars%T_sim-1))
        allocate(sim_result%SWF_sim(pars%T_sim))
        allocate(sim_result%term_sim(pars%T_sim))
        allocate(sim_result%totret_sim(pars%T_sim))

        allocate(sim_result%tau_sim(pars%T_sim,I_par))
        allocate(sim_result%loss_sim(pars%T_sim))
        allocate(sim_result%s_ind_sim(pars%T_sim))
        allocate(sim_result%G_sim(pars%T_sim,I_par))
        allocate(sim_result%G_to_Y_sim(pars%T_sim,I_par))
        allocate(sim_result%Y_sim(pars%T_sim,I_par))
        allocate(sim_result%B_to_GDP_sim(pars%T_sim,I_par))
        allocate(sim_result%B_ex_to_GDP_sim(pars%T_sim))

        !Also initialize these to zero values (so that if I run out of runtime and I want to save at least the
        !series which I managed to compute, I save zeroes and not some nonsense
        sim_result%c_sim = 0.0_dp
        sim_result%l_sim = 0.0_dp
        sim_result%b_sim = 0.0_dp
        sim_result%b_prime_sim = 0.0_dp
        sim_result%a_prime_sim = 0.0_dp
        sim_result%u_c_sim = 0.0_dp
        sim_result%u_l_sim = 0.0_dp
        sim_result%U_sim = 0.0_dp
        sim_result%rho_prime_sim = 0.0_dp
        sim_result%aex_prime_sim = 0.0_dp
        sim_result%b_ex_sim = 0.0_dp
        sim_result%b_ex_prime_sim = 0.0_dp
        sim_result%R_sim = 0.0_dp
        sim_result%SWF_sim = 0.0_dp
        sim_result%term_sim = 0.0_dp
        sim_result%totret_sim = 0.0_dp
        sim_result%tau_sim = 0.0_dp
        sim_result%loss_sim = 0.0_dp
        sim_result%s_ind_sim = 0.0_dp
        sim_result%G_sim = 0.0_dp
        sim_result%G_to_Y_sim = 0.0_dp
        sim_result%Y_sim = 0.0_dp
        sim_result%B_to_GDP_sim = 0.0_dp
        sim_result%B_ex_to_GDP_sim = 0.0_dp
    else if(wtd==0) then
        !Deallocate
        deallocate(sim_result%c_sim)
        deallocate(sim_result%l_sim)
        deallocate(sim_result%b_sim)
        deallocate(sim_result%b_prime_sim)
        deallocate(sim_result%a_prime_sim)
        deallocate(sim_result%u_c_sim)
        deallocate(sim_result%u_l_sim)
        deallocate(sim_result%U_sim)
        deallocate(sim_result%rho_prime_sim)
        deallocate(sim_result%aex_prime_sim)
        deallocate(sim_result%b_ex_sim)
        deallocate(sim_result%b_ex_prime_sim)
        deallocate(sim_result%R_sim)
        deallocate(sim_result%SWF_sim)
        deallocate(sim_result%term_sim)
        deallocate(sim_result%totret_sim)
        deallocate(sim_result%tau_sim)
        deallocate(sim_result%loss_sim)
        deallocate(sim_result%s_ind_sim)
        deallocate(sim_result%G_sim)
        deallocate(sim_result%G_to_Y_sim)
        deallocate(sim_result%Y_sim)
        deallocate(sim_result%B_to_GDP_sim)
        deallocate(sim_result%B_ex_to_GDP_sim)
    else
        write(*,*) 'Error: wrong value of wtd in subroutine sim_result_alloc.'
        sync all
        error stop
    end if

end subroutine sim_result_alloc


!Subroutine sim_series computes simulated series for periods 1,...,T_sim-1
!(so that the total number of simulated elements including thhose for t=0
!will be T_sim)
subroutine sim_series(sim_result,FP_result,pars,grids,x_guess_t1,x_pol_all_unf,sim_shared,s_ind_loaded,sim_ind,CPO_optim_count)
    Use nag_library, Only: e04jyf, nag_wp, e04jcf, e04jcp,e05jaf,e05jbf,E05JBK,e05saf, e05zkf, e05zlf,E05SXM
    !The E05 libraries are not used at this point but may be if we
    !want to look for corner solutions as well (in which case we probably would not have
    !a good initial guess from FP problem and we should use maximization subroutines).

    type(sim_result_type), intent(inout) :: sim_result
    type(FP_result_type), intent(in) :: FP_result
    type(par), intent(in) :: pars
    type(grids_type),intent(in) :: grids
    real(dp), intent(in), dimension(ncvar_par) :: x_guess_t1
    real(dp), dimension(:,:,:,:,:,:), intent(in) :: x_pol_all_unf !this is allocated only if
    !pars%sim_interp = .true.!

    real(dp), dimension(M_par,1) :: ones_col = 1.0_dp

    type(sim_shared_type), codimension[*] :: sim_shared

    integer, intent(in) :: s_ind_loaded(:,:)
    integer, intent(in) :: sim_ind

    integer, dimension(:), allocatable :: CPO_optim_count !a counter of the times when each image found the optimum choice

    integer :: t_ind !period index

    !Pointer to value function. This will be used to make the old value function V_old
    !accessible to the subroutine used in NAG maximization.
    real(dp), dimension(:,:,:,:,:), pointer :: V_old_pntr

    !Values of state variables as saved in the grid (hence the suffix st).
    real(dp), dimension(I_par) :: a_st
    real(dp), dimension(I_par-1) :: rho_st
    real(dp) :: aex_st
    integer :: g_st_ind

    !Variables used in NAG subroutine E04JYF
    real(dp), dimension(ncvar_par) :: xc,xc_last,xc_cent,xc_interp,xc_init,xc_cent_loc !choice variables
    !xc_last contains the optimal choice from the last period. xc_cent is the center
    !around which we perform the permutations.
    !xc interp is the interpolated value
    !It is better to save these separately to avoid confusion
    real(dp) :: fc_last,fc_cent,fc_interp

    real(dp), dimension(ncvar_par) :: rnd_num

    !Best choice so far (used when we are attempting the maximization multiple times)
    real(dp), dimension(ncvar_par) :: bcsf_x !choice
    real(dp) :: bcsf_fc !value (as returned by subroutine, so we are looking for the smallest possible one).
    integer :: att_ind !index for cycling over attempts

    real(dp) :: fc !function value in maximization
    !The following two arrays are arrays of reals and integers which can be used to supply information
    !to the maximized function as an alternative to COMMON block.
    real(dp), dimension(1) :: ruser !Contains no information, only defined for purposes of using a NAG subroutine
    integer, dimension(3) :: iuser !will be used to pass what_binds to some subroutines.
    integer :: ifail, nf
    !Bounds:
    real(dp), dimension(ncvar_par) :: bl
    real(dp), dimension(ncvar_par) :: bu
    real(dp) :: fmax

    integer :: m_ind !for cycling over states

    !Current_period choices (usually state-contingent plans. Only one element of these is saved into sim_series,
    !according to shock realization).
    real(dp), dimension(M_par,I_par) :: u_c,u_l
    real(dp), dimension(M_par,I_par) :: R !I actually compute R for both countries separately to make sure there are no errors.
    !(it should always be equal unless I made an error somewhere).
    real(dp), dimension(M_par,I_par) :: c,l !consumption and labour supply
    real(dp), dimension(M_par,I_par) :: a_pr !next-period adjusted MU-asset holdings
    real(dp), dimension(M_par,1) :: rho_pr !next-period ratio of marginal utilities
    real(dp), dimension(M_par,1) :: aex_pr !next-period MU-adjusted external debt of country 1
    integer :: getavfail

    !Honestly, I am not sure if this makes sense for discarding outliers. If anything, we should use value of
    !loss function (which we can get for every point using function CPC_ret and common block).
    !In either case this is not used for obtaining any of the results which are reported in the paper
    !but for debugging only.
    real(dp), dimension(I_par) :: a_min_slack !slackness of constraint on range of
    !a_prime from below which follows from choice of grid (pars%a_min), for all agents a_prime
    !(positive values mean that the constraint is strictly satisfied, negative mean the opposite, 0 means that the
    !constraint is satisfied with equality).
    real(dp), dimension(I_par) :: a_max_slack !the same but for constraint on a from above.
    real(dp), dimension(:), allocatable :: a_slack_all !contains the lowest slackness of all constraints
    !on agent's a_prime in a period). I.e., the most violated constraint in a period is reflected in this,
    !along with seriousness of the violation.

    integer :: att_counter

    integer, dimension(M_par,4) :: j_guess !The initial guess of position of next-period states
    !in the grid, for evry shock realization.

    integer, dimension(4) :: j_guess_pol !also get j_guess for policy function. In this case only
    !one row - because we do not do one interpolation for every shock realization but one interpolation
    !only (we do not need to compute expectation of policy function as opposed to expectation of value function)

    integer :: attmax_local

    !some other things needed by  e04jyf
    !Try changing the dimension of these things and see what happens (there is only a lower bound)
    integer :: LIW = ncvar_par + 2,LW=max(ncvar_par*(ncvar_par-1)/2+12*ncvar_par,13)
    integer, dimension(ncvar_par + 2) :: IW
    real(dp), dimension(max(ncvar_par*(ncvar_par-1)/2+12*ncvar_par,13)) :: W

    !Common block to be accessible by subroutine CPC_ret. This needs to be exactly the same here
    !and in VFI_solve_gp (and in CPC_ret). Mistakes here are often not detected by the compiler and
    !have catastrophic consequences.

    real(dp) :: deltav_debug,deltax_debug,deltav_debug_avg,deltax_debug_avg !change in optimal policy (the policy
    !that the algorithm finds as optimal after reaching the maximum number of attempts - this is not necessarily the
    !truly optimal policy) over the initial guess expressed as absolute relative deviation. This is saved in deltax.
    !Delta v contains the same information but it is the change in value. The variables with suffix _avg are
    !averages of these over all periods.
    !These are computed and displayed only if pars%sim_debug == .true.
    !value_sum - this is the sum of values of CPC_ret and of current plus continuation value
    !(multiplied by -1 because we are interested in
    real(dp) :: value_sum
    real(dp) :: IG_value !value of initial guess (debug information).
    !The initial guess will be either xc_last, xc_interp, or the better of the two
    !(in case when sim_interp and sim_takebest are both true)
    integer :: PSO_improv_counter !counts the number of times when PSO
    !yielded a better value than the initial guess.
    integer :: takebest_improv_counter !counts the number of times when takebest was used
    !(i.e., when we used policy interpolation but the last-period choice xc_last
    !turned out to be better).

    real(dp), dimension(5,ncvar_par) :: xc_cent_all !this is an array of choices around which we generate the permutations.
    !so far I assume that there are at most 5 different choices around which we permutate. Increase this later if needed.
    integer :: xc_cent_count,xc_cent_ind !the first one is the number of choices that we use. The second one is an index
    !for cycling over these choices.
    logical :: xc_cent_improv2a

    integer :: img_ind !for cycling over images

    !_________The following variables are used only in PSO____________
    integer :: ndim, npar
    integer, parameter :: LIOPTS = 100,lopts = 100
    integer :: iopts(LIOPTS)
    real(dp) :: opts(lopts)
    integer :: itt(6)
    real(dp), dimension(ncvar_par) :: bl_pso
    real(dp), dimension(ncvar_par) :: bu_pso
    integer :: inform
    !__________END of PSO variables____________
    !movea (in the first part of randomization might be a bit risky, careful about it).
    !actually it's always risky. We might come across a point which is a little bit better
    !a few times by chance byt then it might be further from the real optimum than the
    !bounds of the interval - a few unlucky jumps and we mess things up completely! - simulation
    !is very sensitive to small errors.). Move a should always be false, move_b probably as well.
    !these are really risky options!
    logical, parameter :: xc_movea = .false., xc_moveb = .false.

    integer :: bcsf_img !the index of image where we found the best choice (for debug purposes only -
    !for example if this is always 1, we probably made a mistake)

    !Also with PSO - check whether any value of inputs is changes because intent in/out does
    !not work there.
    real(dp) :: cum_SWF !cumulated discounted SWF

    integer :: seed_n
    integer, allocatable :: seed_array(:)


    common /cpc_ret_bloc/ a_st,rho_st,aex_st,g_st_ind,j_guess
    common /V_old_pointer/ V_old_pntr !Block used to pass pointers to value function and grids
    !to some subroutines.



    if(this_image() == 1) then
    !We must be very careful with using variables because the program is run on
    !multiple CPUs, but some variables are only allocated on CPU1 to save memory, and those that are allocated everywhere
    !may have different values between the CPUS.

!_____________________________________________________________________________
    !Save the first-period results computed previously as first elements of the simulated series
    sim_result%c_sim(1,:) = FP_result%c(1,:)
    sim_result%u_c_sim(1,:) = FP_result%u_c(1,:)
    sim_result%u_l_sim(1,:) = FP_result%u_l(1,:)
    sim_result%l_sim(1,:) = FP_result%l(1,:)
    sim_result%b_prime_sim(1,:) = FP_result%b_prime(1,:)
    sim_result%a_prime_sim(1,:) = FP_result%a_prime(1,:)
    sim_result%rho_prime_sim(1) = FP_result%rho_prime
    sim_result%aex_prime_sim(1) = FP_result%aex_prime
    sim_result%b_ex_prime_sim(1) = FP_result%b_ex_prime

    !in the first period I do not save the termination value (we will never
    !do welfare comparisons for one period only anyway).
    sim_result%term_sim = 0.0_dp

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

    !sim_result%s_real_sim(1) = pars%S(1,sim_result%s_ind_sim(1))
    !Draw chain of shock realizations and save it in sim_result
    !(don't save G yet - fill these in gradually).

    !If we load shocks from file, put this into sim_result%s_ind, otherwise generate a random
    !series
    if(pars%load_shocks) then
        sim_result%s_ind_sim(2:pars%T_sim) = s_ind_loaded(2:pars%T_sim,sim_ind)
        !just checking for errors
        if(sim_result%s_ind_sim(1) /= s_ind_loaded(1,sim_ind)) then
            write(*,*) 'Error: in sim_series the first-period shock realization does not agree with the s_ind_loaded.'
            error stop
        end if

    else
    !IF use_P_sim is true, then we do not draw shocks using the transition matrix P, but rather
    !matrix P_sim. This is useful for debugging.
    if(pars%use_P_sim) then
        sim_result%s_ind_sim(2:pars%T_sim) = markov(FP_result%s_ind_init,pars%P_sim,pars%T_sim-1)
    else !this is the default option which should be used all the time outside of debugging
        sim_result%s_ind_sim(2:pars%T_sim) = markov(FP_result%s_ind_init,pars%P,pars%T_sim-1)
    end if
    end if

    !initialize debug variables

    deltax_debug_avg = 0.0_dp
    deltav_debug_avg = 0.0_dp
    value_sum = 0.0_dp
    takebest_improv_counter = 0
    PSO_improv_counter = 0

    xc_cent_all = 0.0_dp

    !In period t=1, the initial guess will be the one obtained by solving the first-period
    !problem for every possible shock realization (to generate a state-contingent plan
    !which is a better guess than a plan which uses the actually realized first period
    !solution in every state). Alternatively (and this is usually better and it is the default
    !option), the initial guess will be obtained using interpolation of the policy
    !function. This is done if pars%interp = .true. This is done later on (and it overwrites
    !the guess set here if that is the case)
    xc_last = x_guess_t1

    end if !(i_img == 1)


    !In periods (indexed from t=2) to T_sim_max, the initial guess will be the optimal solution
    !from the previous period.

    j_guess = 0
    !initialize guess of position of next-period states on the grid. An alternative would be to initialize it
    !again for each max_att (particularly if permutations are large), but it should make very little
    !difference in practice. What is important is that the guess is passed between iterations
    !when maximizing (there can be many iterations, and next-period states should be closely related).
    j_guess_pol = 0


!    write(*,*) 'TEMPORARILY IMPOSING A CUSTOM STATE AT THE BEGINNING OF SIM_SERIES'
!        sim_result%a_prime_sim(1,:) = 0.0_dp
!        sim_result%rho_prime_sim(1) = 1.0_dp
!        sim_result%aex_prime_sim(1) = 0.0_dp
        !sim_result%s_ind_sim(t_ind-1)

    !cycle over time period to generate the simulated series.
    do t_ind=2,pars%T_sim
        xc_cent_count = 0

        if(this_image() == 1) then
        !initialize count of choices. Every period this is set to zero, and if sim_careful,
        !then at several stages a new initial guess is added in xc_cent_all, and the counter is increased.
        xc_cent_improv2a = .false. !this is changed to true if an improvement is found in stage 2a (iterations
        !over sim_attmax2a).

        !At the outset, add xc_last (last period choice) to xc_cent_all.
        if(pars%sim_careful) then
            xc_cent_count = xc_cent_count + 1
            xc_cent_all(xc_cent_count,:) = xc_last
        end if

        !Full debug info
        if(pars%sim_debug) then
            write(*,*) '========Debug in sim_series================='
            write(*,*) 'period t = ',(t_ind-1)
        end if


        !Get states from previous elements of sim_result (so we can use the same code as elsewhere)
        a_st = sim_result%a_prime_sim(t_ind-1,:)
        rho_st = sim_result%rho_prime_sim(t_ind-1)
        aex_st = sim_result%aex_prime_sim(t_ind-1)
        g_st_ind = sim_result%s_ind_sim(t_ind-1)
        !We compute the plan conditional on last-period shock realization s_ind_sim(t-1)!

        !Now image 1 read states from last-period solution. Create a copy of states
        !and pass it to all images.
        sim_shared%a_st = a_st
        sim_shared%rho_st = rho_st
        sim_shared%aex_st = aex_st
        sim_shared%g_st_ind = g_st_ind

        !Initialize best choice found so far. Use the choice from previous period xc_last, or
        !possibly the interpolated policy (this depends on pars%sim_interp, and pars%sim_takebest).

    !if we are using interpolation in simulation, then obtain the initial guess using interpolation
    !of the policy function. Alternatively, use the last-period solution (xc_last).
    if(pars%sim_interp) then

        !So far I do not use an initial guess of position in the grid. This could speed up the interpolation
        !quite a bit possibly as in the case of value function interpolation. I might implement this at some
        !point in the future. For now getting the initial guess using interpolation is not the part
        !where we spend a lot of time - it is the subsequent improvements of this initial guess.

        !if sim_debug, then call a 'debug' version of the interpolating subroutines. These will print
        !a lot of extra information in the terminal.
        if(pars%sim_debug .and. .false.) then !temporarily disable the debug interpolation - too much information.

        !Which interpolation function we use depends on the value of VFI_interp_mode (set in parameter file).
        select case(pars%x_interpolation_mode)
            case(1) !Quadrilinear interpolation
                !This one does not display debug information.
                call interp_x(x_pol_all_unf,[a_st(1),a_st(2),rho_st(1),aex_st],j_guess_pol,min(g_st_ind,pars%M_grid),xc_interp,grids%a1_gr&
                ,grids%a2_gr,grids%rho_gr,grids%aex_gr,size(grids%a1_gr),size(grids%rho_gr),size(grids%aex_gr))
            case(2) !Local Shepard's interpolation
                call interp_x_shep_debug(pars,x_pol_all_unf,[a_st(1),a_st(2),rho_st(1),aex_st],j_guess_pol,min(g_st_ind,pars%M_grid),xc_interp,grids%a1_gr&
                ,grids%a2_gr,grids%rho_gr,grids%aex_gr,size(grids%a1_gr),size(grids%rho_gr),size(grids%aex_gr),&
                pars%shep_n,pars%shep_norm,pars%shep_norm_val,pars%shep_floor)
        end select

        else !do standard interpolation

        !Which interpolation function we use depends on the value of VFI_interp_mode (set in parameter file).
        select case(pars%x_interpolation_mode)
            case(1) !Quadrilinear interpolation
                call interp_x(x_pol_all_unf,[a_st(1),a_st(2),rho_st(1),aex_st],j_guess_pol,min(g_st_ind,pars%M_grid),xc_interp,grids%a1_gr&
                ,grids%a2_gr,grids%rho_gr,grids%aex_gr,size(grids%a1_gr),size(grids%rho_gr),size(grids%aex_gr))
            case(2) !Local Shepard's interpolation
                call interp_x_shep(x_pol_all_unf,[a_st(1),a_st(2),rho_st(1),aex_st],j_guess_pol,min(g_st_ind,pars%M_grid),xc_interp,grids%a1_gr&
                ,grids%a2_gr,grids%rho_gr,grids%aex_gr,size(grids%a1_gr),size(grids%rho_gr),size(grids%aex_gr),&
                pars%shep_n,pars%shep_norm,pars%shep_norm_val,pars%shep_floor)
        end select

        end if !(sim_debug)

        !Now the best choice so far (the only one we tried) and the center for permutations
        !is going to be the interpolated value
        bcsf_x = xc_interp
        call CPC_ret(5,xc_interp,fc_interp,iuser,ruser)
        bcsf_fc = fc_interp
        xc_cent = xc_interp

        !Regardless of whether this is better than the previous best choice - if sim_careful, then
        !save it for generating permutations around it. Even if sim_takebest, we always keep both.
        if(pars%sim_careful) then
            xc_cent_count = xc_cent_count + 1
            xc_cent_all(xc_cent_count,:) = xc_interp
        end if

        !if sim_takebest we compute the values of both xc_last and xc_interp and choose the one
        !with better (smaller) return as bcsf_x, and xc_cent
        if(pars%sim_takebest) then
            call CPC_ret(5,xc_last,fc_last,iuser,ruser)
            if(fc_last < fc_interp) then
                !In the other case the variables already contain those associated with xc_interp
                bcsf_x = xc_last
                bcsf_fc = fc_last
                xc_cent = xc_last
                takebest_improv_counter = takebest_improv_counter + 1
            end if
        end if

    !I need to make sure that I have bcsf_x and xc_center both saved here. no matter what

    else !The case when sim_interp == .false.
        bcsf_x = xc_last
        call CPC_ret(5,bcsf_x,bcsf_fc,iuser,ruser)
        xc_cent = xc_last

        !in this case both the best choice so far and the center for permutations are xc_last
    end if

        !value of initial guess - the value before PSO or any other permutations happen.
        !(the following 2 are used only if sim_debug == .true.)
        IG_value = bcsf_fc
        xc_init = bcsf_x

    !At this stage xc_cent contains the center point for permutations, and bcsf_x contains the
    !best choice so far.

    !Save all the variables which we need to proceed with maximization.
    sim_shared%xc_cent = xc_cent
    sim_shared%xc_cent_all = xc_cent_all
    sim_shared%xc_cent_count = xc_cent_count
    sim_shared%bcsf_x = bcsf_x
    sim_shared%bcsf_fc = bcsf_fc

    end if !(img = 1)


    !Now on all images except 1, read the states and other varaibles needed in maximization.
    !(IMP) - I should perhaps move these repetitive synchronizations into a subroutine
    !(one for syncing states, one for guesses).
    sync all !so other images do not try to read the data before image 1 saved them.
    if(this_image() /= 1) then

     a_st = sim_shared[1]%a_st
     rho_st = sim_shared[1]%rho_st
     aex_st = sim_shared[1]%aex_st
     g_st_ind = sim_shared[1]%g_st_ind

     xc_cent = sim_shared[1]%xc_cent
     xc_cent_all = sim_shared[1]%xc_cent_all
     xc_cent_count = sim_shared[1]%xc_cent_count
     bcsf_x = sim_shared[1]%bcsf_x
     bcsf_fc = sim_shared[1]%bcsf_fc
    end if
    sync all !so image 1 doesn't change these before we are finished reading them


    !Debug: checking input on every image is correct
!    critical
!        write(*,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
!        write(*,*) 'Optimization inputs on image ',this_image()
!        write(*,*) 'a_st = ',a_st
!        write(*,*) 'rho_st = ',rho_st
!        write(*,*) 'aex_st = ',aex_st
!        write(*,*) 'g_st_ind = ',g_st_ind
!
!        write(*,*) 'xc_cent = ',xc_cent
!        write(*,*) 'xc_cent_all = ',xc_cent_all
!        write(*,*) 'xc_cent_count = ',xc_cent_count
!        write(*,*) 'bcsf_x = ',bcsf_x
!        write(*,*) 'bcsf_fc = ',bcsf_fc
!    end critical



        !IF we are only interested in getting a steady state with transfers, skip the maximization.
        !without permutation, and sticking to the previous period choice, this should result
        !in the SS. Otherwise there is an error somewhere.
        if(pars%gen_transfers_ss) then
            if(num_images() > 1) then
                write(*,*) 'gen_transfers_ss when n_images > 1 is not supported. It will probbaly crash. Run this on 1 image.'
            end if
            if(this_image() == 1) then
                goto 545
            else
                return
            end if
        end if

    !Particle swarm optimization.
    if(pars%sim_PSO) then

        ifail = 0
        Call e05zkf('Initialize = E05SAF',iopts,liopts,opts,lopts,ifail)

        !Various parameters may be changed using zkf. See the sample program on website of FL 24.
        !The default values of the following two parameters are E-04. See what different it makes if
        !I set them to E-5 or E-6 in terms of runtime.
        ifail = 0
        Call e05zkf('Local Interior Tolerance = 0.00001',iopts,liopts,opts,lopts,ifail)
        ifail = 0
        Call e05zkf('Local Exterior Tolerance = 0.00001',iopts,liopts,opts,lopts,ifail)

        ndim = 5
        npar = pars%sim_PSO_par !suggested value is 10*ndim (but perhaps I can try more/less).

        ifail = 1 !silent exit (in fact the documentation suggests this value because ifail = 1
        !every time that it is not guaranteed that global minimum was found - which is most of the times).

        !set bounds
        bl_pso = (1-pars%sim_PSO_perm)*xc_cent
        bu_pso = (1+pars%sim_PSO_perm)*xc_cent

        !If any of the bounds are close to zero, do not even attempt PSO.
        !Here I assume that xc_cent is greater than zero (so bl < bu). If not we have other more serious
        !problems, and the algorithm will crash elsewhere anyway.
        if(minval(reshape(abs(bl_pso),[ncvar_par])) < 0.001_dp) then
            goto 111
        end if

        !Call the PSO subroutine
        Call e05saf(ndim,npar,xc,fc,bl_pso,bu_pso,CPC_ret_PSO,E05SXM,iopts,opts, &
            iuser,ruser,itt,inform,ifail)

        !Just in case that the return fc does not correspond to xc (this can sometimes happen if a nag
        !subroutine encounters an error) - recompute the return.
        call CPC_ret(5,xc,fc,iuser,ruser)


        !save the choices into the shared data type
        sim_shared%fc = fc
        sim_shared%xc = xc
        !Now all images finished PSO. sync, find the best choice, and save it as
        !xc (with corresponding return fc) on image 1.

        !Now image 1 cycles over all other images and finds the best choice (minimal)
        sync all
        if(this_image() == 1) then

        do img_ind = 2,num_images()

            !if sim_debug_par, then list the value at every image
            if(pars%sim_debug_par .and. pars%sim_debug) then
                if(img_ind == 2) then !(in the first iter over images first write the info for image 1.
                    write(*,*) '~~~~~~~~par. debug in simulation (PSO stage)~~~~~~~~~'
                    write(*,*) 'Listing optimal value and choice found by every image...'
                    write(*,*) 'i_img = ',1,'fc = ',fc,'x = ',xc
                    write(*,*) ' ' !empty row (kinda needed due to xc, this will make it messy).
                end if
                    write(*,*) 'i_img = ',img_ind,'fc = ',sim_shared[img_ind]%fc,'x = ',sim_shared[img_ind]%xc
                    write(*,*) ' ' !empty row
            end if

            if(sim_shared[img_ind]%fc < fc) then
                fc = sim_shared[img_ind]%fc
                xc = sim_shared[img_ind]%xc
            end if
        end do
        end if
        sync all !make sure data on other images is not changed before we read it.



        !This part determines what we do with the choice found (basically xc_cent_all and xc_cent determination,
        !and it also compute some debug info).
        !only on image 1


        if(this_image() == 1) then

        !now if this is better (fc < bcsf_fc) save xc as bcsf_x and as xc_cent (so that permutations are done around
        !this point). This could actually potentially lead to some instability.
        if(fc < bcsf_fc) then
            if(pars%sim_debug) then
                write(*,*) 'PSO improvement over init. guess in this iteration = ',bcsf_fc - fc
                PSO_improv_counter = PSO_improv_counter + 1
            end if
            bcsf_x = xc
            bcsf_fc = fc
            xc_cent = xc
        end if

        !If sim_careful, save this for later use regardless of whether it is better than
        !the best choice so far
        if(pars%sim_careful) then
            xc_cent_count = xc_cent_count + 1
            xc_cent_all(xc_cent_count,:) = xc
        end if

        !If PSO causes problems then manually check that it did not change any of its inputs
        !which are meant for use with other subroutines further below.
        !It should not but based on past experience it could feasibly happen.
        sim_shared%xc_cent = xc_cent
        sim_shared%xc_cent_all = xc_cent_all
        sim_shared%xc_cent_count = xc_cent_count
        sim_shared%bcsf_x = bcsf_x
        sim_shared%bcsf_fc = bcsf_fc
        end if !(i_img == 1)

        !Now we need to syncrhonize the guesses and bcsf again.

    sync all !so other images do not try to read the data before image 1 saved them.
    if(this_image() /= 1) then
     xc_cent = sim_shared[1]%xc_cent
     xc_cent_all = sim_shared[1]%xc_cent_all
     xc_cent_count = sim_shared[1]%xc_cent_count
     bcsf_x = sim_shared[1]%bcsf_x
     bcsf_fc = sim_shared[1]%bcsf_fc
    end if
    sync all !so image 1 doesn't change these before we are finished reading them
    !(only do this synchronization if PSO was used - otherwise there was no change since the last sync)

    111 end if !end of PSO

        !this should not be necessary (but no harm done). Test removing it later.
        call CPC_ret(5,bcsf_x,bcsf_fc,iuser,ruser)

        !Now the first part of search. Take pars%attmax2a permutations around xc_center.
        !evaluate each point and if it's better than the previous choice, save it.
        !Also if pars%xc_center_move, then use that point as the new center for optimization.
        if(t_ind == 2) then
            attmax_local = 2*pars%sim_attmax2a
        else
            attmax_local = pars%sim_attmax2a
        end if

        !If sim_par_keepattmax, then we do not change the number of attempts depending
        !on the number of CPUs. Otherwise change it so as to keep the total number of attempts
        !approximately unchanged.
        if(.not. pars%sim_par_keepattmax) then
            if(real(attmax_local,dp)/real(num_images(),dp) > 1.0E-15_dp) then
            !(treat the zero case separately due to possible rounding errors)
                attmax_local = ceiling(real(attmax_local,dp)/real(num_images(),dp))
            else
                attmax_local = 0
            end if
        end if



        xc_cent_ind = 1
        do att_ind = 1,attmax_local
                !If sim_careful, alternate between all xc_cent_all permutation centers. Othwerwise
                !just use xc_cent
                if(pars%sim_careful) then
                    xc_cent_loc = xc_cent_all(xc_cent_ind,:)
                    !if we reached the end then reset the counter. Otherwise increment the counter by 1
                    if(xc_cent_ind == xc_cent_count) then
                        xc_cent_ind = 1
                    else
                        xc_cent_ind = xc_cent_ind + 1
                    end if
                else
                    xc_cent_loc = xc_cent
                end if


                !Draw a random point close to xc_cent (permutations driven by sim_perm2).
                !(it's a box, uniform draw).
                call random_number(rnd_num)
                rnd_num = rnd_num * pars%sim_perm2 * 2.0_dp + 1.0_dp - pars%sim_perm2
                xc = xc_cent_loc * rnd_num

                !Evaluate the point
                call CPC_ret(5,xc,fc,iuser,ruser)

                !If it is better than the best one found previously, save it.
                !If xc_movea, then also make xc_cent the new center of permutation,
                !otherwise it will be moved only once at the very end (before the next stage).
            if(fc < bcsf_fc) then
                bcsf_x = xc
                bcsf_fc = fc
                if(xc_movea) xc_cent = xc
                xc_cent_improv2a = .true.
            end if

        end do

        !Now every image saves its best choice found so far
        sim_shared%bcsf_x = bcsf_x
        sim_shared%bcsf_fc = bcsf_fc
        sync all
        !Image 1 reads the data and figures out whether some of the choice was better
        if(this_image() == 1) then
        do img_ind = 2,num_images()

            !if sim_debug_par, then list the value at every image
            if(pars%sim_debug_par .and. pars%sim_debug) then
                if(img_ind == 2) then !(in the first iter over images first write the info for image 1.
                    write(*,*) '~~~~~~~~par. debug in simulation (stage 2a)~~~~~~~~~'
                    write(*,*) 'Listing optimal value and choice found by every image...'
                    write(*,*) 'i_img = ',1,'fc = ',bcsf_fc,'x = ',bcsf_x
                    write(*,*) ' ' !empty row (kinda needed due to xc, this will make it messy).
                end if
                    write(*,*) 'i_img = ',img_ind,'fc = ',sim_shared[img_ind]%bcsf_fc,'x = ',sim_shared[img_ind]%bcsf_x
                    write(*,*) ' ' !empty row
            end if


            if(sim_shared[img_ind]%bcsf_fc < bcsf_fc) then
                bcsf_fc = sim_shared[img_ind]%bcsf_fc
                bcsf_x = sim_shared[img_ind]%bcsf_x
                xc_cent_improv2a = .true.
            end if
        end do

            !Also if sim_debug display a message that stage 2a found an improvement.

        !the best choice found so far will be the new center of permutations. This is actually
        !quite problematic I think if the number of attempts is not great. But if sim_careful is true,
        !then it should not matter.
        xc_cent = bcsf_x

        !Handle this carefully. Basicaly on image 1 we know whether improvement happened (xc_cent_improv2a = .true.)
        !when reading data from other images, set xc_cent_improv2a to true if bcsf from that image was better.
        !If we are 'careful' then we add the old best choice to the list of permutation centers (but only
        !if some improvement actually happened).
        if(pars%sim_careful .and. xc_cent_improv2a) then
            xc_cent_count = xc_cent_count + 1
            xc_cent_all(xc_cent_count,:) = bcsf_x
        end if
        !img 1 saved all the maximization data for sending to other images:
        sim_shared%xc_cent = xc_cent
        sim_shared%xc_cent_all = xc_cent_all
        sim_shared%xc_cent_count = xc_cent_count
        sim_shared%bcsf_x = bcsf_x
        sim_shared%bcsf_fc = bcsf_fc
        end if !(img = 1)

    !images 2,...,N read data from image 1
    sync all !so other images do not try to read the data before image 1 saved them.
    if(this_image() /= 1) then
     xc_cent = sim_shared[1]%xc_cent
     xc_cent_all = sim_shared[1]%xc_cent_all
     xc_cent_count = sim_shared[1]%xc_cent_count
     bcsf_x = sim_shared[1]%bcsf_x
     bcsf_fc = sim_shared[1]%bcsf_fc
    end if
    sync all !so image 1 doesn't change these before we are finished reading them
    !(only do this synchronization if PSO was used - otherwise there was no change since the last sync)



        !In the second period (t=1) the previous solution might not be as good a guess
        !because it was obtained using a different optimization subroutine, etc. So use twice as
        !many attempts in the second period to account for this.
        if(t_ind == 2) then
            attmax_local = 2*pars%sim_attmax2b
        else
            attmax_local = pars%sim_attmax2b
        end if
        if(.not. pars%sim_par_keepattmax) then
            if(real(attmax_local,dp)/real(num_images(),dp) > 1.0E-15_dp) then
            !(treat the zero case separately due to possible rounding errors)
                attmax_local = ceiling(real(attmax_local,dp)/real(num_images(),dp))
            else
                attmax_local = 0
            end if
        end if

        xc_cent_ind = 1
        !Attempt the maximization pars%sim_attmax2 times.
        do att_ind = 1,attmax_local
                if(pars%sim_careful) then
                    xc_cent_loc = xc_cent_all(xc_cent_ind,:)
                    !if we reached the end then reset the counter. Otherwise increment the counter by 1
                    if(xc_cent_ind == xc_cent_count) then
                        xc_cent_ind = 1
                    else
                        xc_cent_ind = xc_cent_ind + 1
                    end if
                else
                    xc_cent_loc = xc_cent
                end if


            !Generate a permutation of the initial guess around xc_last (previous-period choice).
            !Only do this if this is not the first attempt, otherwise use the previous choice directly.
            !Also if we are careful, the take an unpermutated maximization attempt at each of the xc_cent_count points.
            if( ((att_ind > 1 .and. .not. pars%sim_careful) .or. &
            (att_ind > xc_cent_count .and. pars%sim_careful)) .or. this_image() /= 1 ) then
                !rnd_num contains 5 numbers to be used for permutation. of xc_init.
                call random_number(rnd_num)
                rnd_num = rnd_num * pars%sim_perm2 * 2.0_dp + 1.0_dp - pars%sim_perm2
                xc = xc_cent_loc * rnd_num
            else
                xc = xc_cent_loc
            end if

            !Perform the maximization, starting at initial guess xc (which is the permutated one)
            ifail = 1 !silent exit
            bl = 0.0001_dp !Set lower bounds
            bu = pars%x_bu
            select case(pars%opt_subroutine)
                case(1) !Quadratic approximation method
                    call e04jcf(CPC_ret2,ncvar_par,ncvar_par*2+1,xc,bl,bu,0.05_dp,&
                    pars%sim_rhoend,e04jcp,pars%sim_maxiter,fc,nf,iuser,ruser,ifail)
                case(2) !Newton-like method - this tends to lead to crashes on HPC for some reason.
                    !The parameters are not quite correct (because the number of choice variables is now different).
                    call e04jyf(ncvar_par,0,CPC_ret,bl,bu,xc,fc,iw,liw,w,lw,iuser,ruser,ifail)
                case default
                    write(*,*) 'Error in subroutine sim_series - wrong value of optim_subroutine'
                    error stop
            end select

            !If the choice is better than the previous best one, save it
            if(fc < bcsf_fc) then
                bcsf_fc = fc
                bcsf_x = xc
                !If move_b then also use it as the new center for permutations
                if(xc_moveb) xc_cent = xc
            end if
        end do
        !Update the last period optimal choice (to get initial guess the next period).

        !Now every image should have the best solution (for this period) in bcsf_fc, bcsf_x,
        !so copy it to shared data, sync, and read it on image 1.
        sim_shared%bcsf_x = bcsf_x
        sim_shared%bcsf_fc = bcsf_fc
        sync all

        !Image 1 reads the data and figures out whether some of the choice was better
        if(this_image() == 1) then
        bcsf_img = 1
        do img_ind = 2,num_images()

            !if sim_debug_par, then list the value at every image
            if(pars%sim_debug_par .and. pars%sim_debug) then
                if(img_ind == 2) then !(in the first iter over images first write the info for image 1.
                    write(*,*) '~~~~~~~~par. debug in simulation (stage 2b)~~~~~~~~~'
                    write(*,*) 'Listing optimal value and choice found by every image...'
                    write(*,*) 'i_img = ',1,'fc = ',bcsf_fc,'x = ',bcsf_x
                    write(*,*) ' ' !empty row (kinda needed due to xc, this will make it messy).
                end if
                    write(*,*) 'i_img = ',img_ind,'fc = ',sim_shared[img_ind]%bcsf_fc,'x = ',sim_shared[img_ind]%bcsf_x
                    write(*,*) ' ' !empty row
            end if

            if(sim_shared[img_ind]%bcsf_fc < bcsf_fc) then
                bcsf_fc = sim_shared[img_ind]%bcsf_fc
                bcsf_x = sim_shared[img_ind]%bcsf_x
                bcsf_img = img_ind
            end if
        end do

        if(pars%sim_debug) then
        write(*,*) '_____________________________________'
        write(*,*) 'States at the beginning of the period:'
        write(*,*) 'a_st = ',a_st
        write(*,*) 'aex_st = ',aex_st
        write(*,*) 'rho_st = ',rho_st
        write(*,*) 'g_st_ind = ',g_st_ind

        write(*,*) '_____________________________________'
        write(*,*) 'Initial guess:'
        write(*,*) 'xc = ',xc_last
        write(*,*) 'value (-1) = ',(IG_value*(-1.0_dp))
        write(*,*) '_____________________________________'

        write(*,*) 'Optimal choice:'
        write(*,*) 'xc = ',bcsf_x


        deltax_debug = sum(reshape(abs(bcsf_x - xc_last)/(abs(xc_last) + 1.0_dp),[ncvar_par]))/real(ncvar_par,dp)
        write(*,*) 'delta x (avg. abs. rel. dev.) = ', deltax_debug
        !(averaged across the 5 choice variables, the absolute relative deviation).
        write(*,*) 'value (-1) = ',(bcsf_fc*(-1.0_dp))
        deltav_debug = (bcsf_fc*(-1.0_dp) - IG_value*(-1.0_dp))

        write(*,*) 'Increase of value over I.G. = ',deltav_debug

        !Also contribute to the averages.
        deltax_debug_avg = deltax_debug_avg + deltax_debug/real(pars%T_sim-1,dp) !T-1 because we do not compute this in the first period
        deltav_debug_avg = deltav_debug_avg + deltav_debug/real(pars%T_sim-1,dp)
        !But I also need to know the sum of value (that is the best indicator of performance_
        value_sum = value_sum - bcsf_fc

        write(*,*) 'Optimal choice found by image ',bcsf_img

        write(*,*) '============================================'
        !end of debug
        end if

        !increase count of optimal points found for the image in question
        if(pars%sim_debug_par) then
            CPO_optim_count(bcsf_img) = CPO_optim_count(bcsf_img) + 1
        end if

        xc_last = bcsf_x

        !now bscf_x contains the best choice found (in vector form). We need to recover
        !all the remaining variables using subroutine get_av, and save these into the appropriate
        !elements of sim_series
        545 call get_av(pars,bcsf_x,a_st,rho_st,aex_st,g_st_ind,getavfail,c,l,a_pr,rho_pr,aex_pr,.false.)

        !Also get marginal utilities of consumption and labour in both countries. This could eventually be added to get_av
        !but usually when using the subroutne, we do not need these.
        call util_c(c,u_c,pars%A_par,pars%sigma_par)
        call util_l(l,u_l,pars%B_par,pars%Gamma_par)

        !We got a state-contingent plan for all the variables. We only save the row corresponding to the
        !actualy realized state (which is sim_result%s_ind_sim(t_ind)).

        !We've got plans for most of the variables which we save. The current-period
        !shock index is sim_series%s_ind_sim(t_ind). Take that row of the plans and save it in
        !sim_series (so we save the actually realized values)

        !First save all the variables which will be states the next period (we already generated
        !shock series at the outset)
        sim_result%a_prime_sim(t_ind,:) = a_pr(sim_result%s_ind_sim(t_ind),:)
        sim_result%rho_prime_sim(t_ind) = rho_pr(sim_result%s_ind_sim(t_ind),1)
        sim_result%aex_prime_sim(t_ind) = aex_pr(sim_result%s_ind_sim(t_ind),1)

        !We have saved the variables essential for solving the model (states)
        !Now we can compute other things and save them into sim_result
        !(government expenditure in levels, taxes, etc.)
        sim_result%c_sim(t_ind,:) = c(sim_result%s_ind_sim(t_ind),:)
        sim_result%l_sim(t_ind,:) = l(sim_result%s_ind_sim(t_ind),:)
        sim_result%u_c_sim(t_ind,:) = u_c(sim_result%s_ind_sim(t_ind),:)
        !Asset holdings (without MU adjustment) - subroutine get_av doesn't return these, we need to compute them
        !manually
        sim_result%b_prime_sim(t_ind,:) = pars%beta * sim_result%a_prime_sim(t_ind,:) / sim_result%u_c_sim(t_ind,:)
        sim_result%b_ex_prime_sim(t_ind) = pars%beta * sim_result%aex_prime_sim(t_ind) / sim_result%u_c_sim(t_ind,1)
        !Also save b_ex_prime.
        sim_result%u_l_sim(t_ind,:) = u_l(sim_result%s_ind_sim(t_ind),:)

        !Get (gross) rate of return on the bonds
        !in periods t, when we know period t consumption plan, we compute period
        !t-1 rate of return. The length of the simulated series for R_t will be
        !one less than the length of the other series. Conceptually we could easily
        !get the missing element but it is cumbersome.

        !Just checking - if the interest rate in the countries is not the same, write a warning about this,
        !and display the interest rates. This should never happen because even if the choice (policy)
        !is suboptimal, interest rates should be the same due to the rho = Eu1/Eu2 constraint
        !Interest rate computed using MU_C of country 1
        R(1,1) = (1.0_dp/pars%beta)*sim_result%u_c_sim(t_ind-1,1)/&
            maxval(matmul(pars%P(sim_result%s_ind_sim(t_ind-1):sim_result%s_ind_sim(t_ind-1),:),u_c(:,1:1)))

        !Interest rate computed using MU_C of country 2
        R(1,2) = (1.0_dp/pars%beta)*sim_result%u_c_sim(t_ind-1,2)/&
            maxval(matmul(pars%P(sim_result%s_ind_sim(t_ind-1):sim_result%s_ind_sim(t_ind-1),:),u_c(:,2:2)))

        if(abs(R(1,1) - R(1,2)) > 0.0000000000001_dp) then
            write(*,*) '______________________________________'
            write(*,*) 'Warning: interest rate is not always equal between the countries.'
            write(*,*) 'Difference = ',abs(R(1,1) - R(1,2))
            write(*,*) '______________________________________'
        end if

        !Average of the interest rates in the two countries. These should be exactly the same,
        !in practice they differ by about 1E-16 due to numerical errors.
        sim_result%R_sim(t_ind-1) = (R(1,1) + R(1,2)) / 2.0_dp

        !Compute the terminal value - this is the actually realized continuation value in this
        !period, discounted to this period. This is useful in welfare analysis if we want to
        !cut off at a particular date.
        !In the case of SS with transfers this is compute manually because we do not have the value function
        !in this case.

        if(.not. pars%gen_transfers_ss) then
            !the value of iuser tells the program what to do (see description of subroutine CPC_ret_adv)
            !We want the value function for the last period.
            iuser = [2,sim_result%s_ind_sim(t_ind),0]
            call CPC_ret_adv(ncvar_par,xc_last,sim_result%term_sim(t_ind),iuser,ruser)
            !The subroutine returns minus 1 times the continuation value...
            sim_result%term_sim(t_ind) = -sim_result%term_sim(t_ind)
            iuser = 0 !reset iuser back to avoid errors if it is used again
        end if

        end if !img = 1
    end do

    !End of the essential part of simulation. All we do now is computing some more statistics.
    !(and we do it on image 1 only)
    if(this_image() /=1) return

    !display more debug info 
    if(pars%sim_debug) then
        write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        write(*,*) '(The averages are over time excluding the first period problem)'
        write(*,*) 'AVERAGE delta x avg = ', deltax_debug_avg
        write(*,*) 'AVERAGE Increase of value over I.G.= ',deltav_debug_avg

        write(*,*) 'Total value (-CPC_ret sum) = ',value_sum
        write(*,*) 'Average Total value = ',(value_sum/real(pars%T_sim-1,dp))

        write(*,*) 'Share where PSO improv. occured = ',real(pso_improv_counter,dp)/real(pars%T_sim-1,dp)
        write(*,*) 'Share where xc_last better than xc_interp = ',real(takebest_improv_counter,dp)/real(pars%T_sim-1,dp)
        write(*,*) '(note : these are computed only if sim_PSO and sim_takebest are .true.).'
        write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    end if

    !Now we have some series saved in sim_result. We can compute a few others
    sim_result%b_sim(1,:) = pars%b_init(1,:) !assets including the initial holdings
    sim_result%b_sim(2:pars%T_sim,:) = sim_result%b_prime_sim(1:pars%T_sim-1,:)

    sim_result%b_ex_sim(1) = pars%b_ex_init !assets including the initial holdings
    sim_result%b_ex_sim(2:pars%T_sim) = sim_result%b_ex_prime_sim(1:pars%T_sim-1)

    cum_SWF = 0.0_dp !initite the cumulative SWF

    do t_ind=1,pars%T_sim
        !realized government expenditure in all countries

        !If we are getting the SS with transfers, use G_copy to get the government expenditures
        !(because the matrix G was overridden at the outset to generate a simulation
        !corresponding to this SS)
        if(pars%gen_transfers_ss) then
            sim_result%G_sim(t_ind,:) = pars%G_copy(sim_result%s_ind_sim(t_ind),:)
        else !the default option (almost always the case)
            sim_result%G_sim(t_ind,:) = pars%G(sim_result%s_ind_sim(t_ind),:)
        end if

        !Get the tax rates in all countries
        sim_result%tau_sim(t_ind,:) = sim_result%u_l_sim(t_ind,:)/(sim_result%u_c_sim(t_ind,:)&
        *pars%theta_0(1,:)) + 1.0_dp

        !output in both countries (aggregate, not per-capita). So multiply by mass
        sim_result%Y_sim(t_ind,:) = pars%mass(1,:) * sim_result%l_sim(t_ind,:) * pars%theta_0(1,:)

        !government debt to GDP
        sim_result%B_to_GDP_sim(t_ind,:) = sim_result%b_sim(t_ind,:)/sim_result%Y_sim(t_ind,:)

        !Net external debt to GDP (of country 1)
        sim_result%B_ex_to_GDP_sim(t_ind) = sim_result%b_ex_sim(t_ind)/sim_result%Y_sim(t_ind,1)

        !Utility in both countries (the actual realized utility, not expected utility)
        call util(sim_result%c_sim(t_ind:t_ind,:),sim_result%l_sim(t_ind:t_ind,:),sim_result%U_sim(t_ind:t_ind,:),&
        pars%A_par,pars%B_par,pars%sigma_par,pars%Gamma_par)

        !Now get value of SWF (the actual realized one, not expected)
        sim_result%SWF_sim(t_ind:t_ind) = reshape(matmul(sim_result%U_sim(t_ind:t_ind,:) &
            * pars%mass * pars%alpha,ones_col_I),[1])

        !Now, if we are generating steady state with transfers, set the termination
        !return to sum of series of SWF in period 1 consumed forever, discounted by beta.
        !Do this only once because it is constant over time.
        if(pars%gen_transfers_ss .and. t_ind == 1) then
            sim_result%term_sim = (pars%beta/(1.0_dp - pars%beta)) * sim_result%SWF_sim(t_ind)
        end if

        !And get the total return up unto the period (if we stopped right now). This is useful
        !mainly for debugging simulations and ultimately comparing how well two different solutions
        !of the same problem fare. In the fisrt period it does not contain the continuation return.
        !add the actually realized SWF(discounted)
        cum_SWF = cum_SWF + sim_result%SWF_sim(t_ind) * (pars%beta**(t_ind -1))

        sim_result%totret_sim(t_ind) = cum_SWF + sim_result%term_sim(t_ind) * (pars%beta**(t_ind -1))

    end do

    !The following is for discarding 'outliers' = bad solutions which diverge. This is to be used with caution
    !if disc_share = 0 in parameter file, no discarding happens.        
    
    !It is also useful for detecting errors.
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

    !Now compute the average a_slack_all. We could truncate these from above by 0, so that
    !violations of constraints would not be countered by strict satisfaction of the constraints
    !in later periods. However, we're using this to eliminate outliers, and not truncating this will
    !favour solutions not only violating constraints but also staying further from the bounds most
    !of the time (other things being equal), so it's probably better to leave it like this.
    !But then this value can't be interpreted in a straightforward way and just because it's not
    !negative doesn't mean things are fine. Instead use quantiles of a_prime for diagnostics (which are also saved).
    sim_result%a_slack_all_avg = 0.0_dp
    do t_ind = 1,pars%T_sim
        sim_result%a_slack_all_avg = sim_result%a_slack_all_avg + a_slack_all(t_ind)/real(pars%T_sim,dp)
    end do

    deallocate(a_slack_all)

end subroutine sim_series


!Subroutine sim_stat computes various statistics from multiple simulations (computed beforehand).
subroutine sim_stat(sim_res_all,pars,sim_avg,N_sim_done)
    use nag_library, Only: G01AMF

    type(sim_result_all_type), intent(inout) :: sim_res_all !on the input contains all simulated series,
    !some additional results will be added (such as residuals, statistics, etc.)
    type(par), intent(in) :: pars
    type(sim_result_type), intent(inout) :: sim_avg !will contain the sample average of all series
    integer, intent(in) :: N_sim_done

    integer :: sim_ind

    type(par) :: pars_comm !copy of pars to be passed to subroutines for manipulating series
    !which are called using operator overloading (hence can't be passed as an input directly
    !or at least I don't know how to do that)

    integer :: m_ind,n_ind,t_ind !indices for cycling over shock realizations, samples, and periods

    integer :: num_of_m !number of realizations of a given shock index (across all N_sim*T_sim observations).

    real(dp), allocatable, dimension(:) :: tmp_vector !for computing quantiles

    !some variables used by NAG subroutine used to compute quantiles:
    integer :: n,nq,ifail
    real(kind = nag_wp) :: Q(2),QV(2)

    !This is used just so we can use operator overloading
    common /pars_comm/ pars_comm

    pars_comm = pars

    !First compute the sample average series.
    !(this uses operator overloading)
    sim_avg = sim_div(sim_res_all%SRA(1),real(N_sim_done,dp))

    do sim_ind = 2,N_sim_done
        sim_avg = sim_avg + sim_div(sim_res_all%SRA(sim_ind),real(N_sim_done,dp))
    end do

    !Compute average aggregate product in all countries for each shock realization (and set it to 0.0 for shock
    !realizations which are never attained). This is useful for calibration (particularly when shocks
    !are aggregate).
    do m_ind = 1,M_par
        num_of_m = 0
        sim_res_all%avg_y(m_ind,:) = 0.0_dp !initialize

        !cycle over all samples and all periods and add the product to sim_res%y_avg(m_ind)
        do n_ind = 1,N_sim_done
            do t_ind = 1,pars%T_sim
                !if the shock realizations for this observation has index m_ind, add it to the
                !average.
                if(sim_res_all%SRA(n_ind)%s_ind_sim(t_ind) == m_ind) then
                    sim_res_all%avg_y(m_ind,:) = sim_res_all%avg_y(m_ind,:) + sim_res_all%SRA(n_ind)%Y_sim(t_ind,:)
                    num_of_m = num_of_m + 1
                end if
            end do
        end do
    !divide by the number of observations to get the sample mean
    if(num_of_m > 0) then
        sim_res_all%avg_y(m_ind,:) = sim_res_all%avg_y(m_ind,:)/real(num_of_m,dp)
    else
        sim_res_all%avg_y(m_ind,:) = 0.0_dp !Just in case that some shock realization is never attained.
    end if
    end do

    !Also compute quantiles for state variables - for debugging purposes (detection of cases
    !where bounds are often hit or closely approached for some shock realizations, even though
    !the average is not close to the bounds).

    !The quantiles are computed for every period (across samples)
    allocate(tmp_vector(pars%N_sim))

    !get the quantiles for a_prime in country 1
    do t_ind = 1,pars%T_sim
        !cycle over samples and save all values of a_prime into vector
        do n_ind = 1,N_sim_done
            tmp_vector(n_ind) = sim_res_all%SRA(n_ind)%a_prime_sim(t_ind,1)
        end do

        N = N_sim_done !number of points in the vector
        nq = 2 !number of quantiles
        q = [0.05_dp, 0.95_dp] !quantile values
        qv = 0.0_dp
        ifail = -1

        call G01AMF(N,tmp_vector,NQ,Q,QV,IFAIL)

        sim_res_all%a_prime_05(t_ind,1) = qv(1)
        sim_res_all%a_prime_95(t_ind,1) = qv(2)
    end do

    !get the quantiles for a_prime in country 2
    do t_ind = 1,pars%T_sim
        !cycle over samples and save all values of a_prime into vector
        do n_ind = 1,N_sim_done
            tmp_vector(n_ind) = sim_res_all%SRA(n_ind)%a_prime_sim(t_ind,2)
        end do

        N = N_sim_done !number of points in the vector
        nq = 2 !number of quantiles
        q = [0.05_dp, 0.95_dp] !quantile values
        qv = 0.0_dp
        ifail = -1

        call G01AMF(N,tmp_vector,NQ,Q,QV,IFAIL)

        sim_res_all%a_prime_05(t_ind,2) = qv(1)
        sim_res_all%a_prime_95(t_ind,2) = qv(2)
    end do

    !IF there are more than 2 countries, display warning
    if(I_par > 2) then
        write(*,*) 'Warning: Subroutine sim_stat needs to be generalized for I_par > 2.'
        write(*,*) '(program execution will not be stopped but some results will be nonsensical)'
    end if

    !get the quantiles for rho_prime
    do t_ind = 1,pars%T_sim
        !cycle over samples and save all values of a_prime into vector
        do n_ind = 1,N_sim_done
            tmp_vector(n_ind) = sim_res_all%SRA(n_ind)%rho_prime_sim(t_ind)
        end do

        N = N_sim_done !number of points in the vector
        nq = 2 !number of quantiles
        q = [0.05_dp, 0.95_dp] !quantile values
        qv = 0.0_dp
        ifail = -1

        call G01AMF(N,tmp_vector,NQ,Q,QV,IFAIL)

        sim_res_all%rho_prime_05(t_ind) = qv(1)
        sim_res_all%rho_prime_95(t_ind) = qv(2)
    end do

     !get the quantiles for aex_prime
    do t_ind = 1,pars%T_sim
        !cycle over samples and save all values of a_prime into vector
        do n_ind = 1,N_sim_done
            tmp_vector(n_ind) = sim_res_all%SRA(n_ind)%aex_prime_sim(t_ind)
        end do

        N = N_sim_done !number of points in the vector
        nq = 2 !number of quantiles
        q = [0.05_dp, 0.95_dp] !quantile values
        qv = 0.0_dp
        ifail = -1

        call G01AMF(N,tmp_vector,NQ,Q,QV,IFAIL)

        sim_res_all%aex_prime_05(t_ind) = qv(1)
        sim_res_all%aex_prime_95(t_ind) = qv(2)
    end do

    deallocate(tmp_vector)

end subroutine sim_stat



subroutine save_series(sim_result,folder_name,pars,N_sim_done)
    use ifport !so we can use system()

    type(sim_result_type) :: sim_result
    character*(*), intent(in) :: folder_name
    type(par), intent(in) :: pars
    integer, intent(in) :: N_sim_done

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

    file_path = 'results/'//trim(folder_name)//'/sim/s_ind.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%s_ind_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/R.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%R_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/tau.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%tau_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/g.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%g_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/rho_prime.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%rho_prime_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/aex_prime.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%aex_prime_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/Y.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%Y_sim
    close(unit = 20)


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

    file_path = 'results/'//trim(folder_name)//'/sim/B_to_GDP.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%B_to_GDP_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/b_ex.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%b_ex_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/B_ex_to_GDP.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%B_ex_to_GDP_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/u.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%u_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/SWF.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%SWF_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/term.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%term_sim
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/totret.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_result%totret_sim
    close(unit = 20)

    !(3)Also save M_par and I_par and T_sim, N_sim (for Matlab to load) in one vector : matlab_info.out
    file_path = 'results/'//trim(folder_name)//'/sim/matlab_control.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) M_par
    write(unit=20, fmt=*) I_par
    write(unit=20, fmt=*) pars%T_sim
    write(unit=20, fmt=*) N_sim_done !save the actual number of simulations that were completed in time
    close(unit = 20)

    !If we got here without crash it means that the results were saved successfuly
    write(*,*) 'Simulated series were saved in folder results/',trim(folder_name),'.'

end subroutine save_series

!Subroutine save_stat saves some statistics about the solution (which are saved in
!type sim_series_all, along with all the simulated series).
subroutine save_stat(sim_res_all,folder_name,pars,N_sim_done)
    type(sim_result_all_type), intent(in) :: sim_res_all
    character*(*), intent(in) :: folder_name
    type(par), intent(in) :: pars
    integer, intent(in) :: N_sim_done

    character*256 :: shell_cmd,file_path
    integer(4) :: sys_call_result
    logical :: folder_exists

    integer :: sim_ind

    real(dp), dimension(N_sim_done,pars%T_sim,I_par) :: tmp_array_i !for saving all series which have 1 element per country
    real(dp), dimension(N_sim_done,pars%T_sim) :: tmp_array !for saving all series which have 1 element


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

    file_path = 'results/'//trim(folder_name)//'/sim/stat_avg_y.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%avg_y
    close(unit = 20)

    !Also save Pareto weight on agent 1 ('rich') - this is useful for plotting
    !various results as a function of the paramter.
    !(it could be read directly from parameter file but it's easier this way)
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

    file_path = 'results/'//trim(folder_name)//'/sim/stat_aex_prime_05.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%aex_prime_05
    close(unit = 20)

    file_path = 'results/'//trim(folder_name)//'/sim/stat_aex_prime_95.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) sim_res_all%aex_prime_95
    close(unit = 20)


    !For most variables, we want to save all simulated series (not just the mean)
    !This is useful for generating fan charts (which are both interesting to put in the paper and useful
    !for debugging), and for computing various statistics in Matlab (more convenient than doing this in Fortran).

    !All these series are saves with the same file name as the average series, with the prefix 'all_'

    !save the tax rate for all series
    do sim_ind = 1,N_sim_done
        tmp_array_i(sim_ind,:,:) = sim_res_all%SRA(sim_ind)%tau_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/all_tau.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array_i
    close(unit = 20)

    !save the government debt for all series
    do sim_ind = 1,N_sim_done
        tmp_array_i(sim_ind,:,:) = sim_res_all%SRA(sim_ind)%b_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/all_b.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array_i
    close(unit = 20)

    !save the external debt of country 1 for all series (1 dim series)
    do sim_ind = 1,N_sim_done
        tmp_array(sim_ind,:) = sim_res_all%SRA(sim_ind)%b_ex_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/all_b_ex.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array
    close(unit = 20)

    !Also save debt and external debt relative to GDP
    do sim_ind = 1,N_sim_done
        tmp_array_i(sim_ind,:,:) = sim_res_all%SRA(sim_ind)%b_to_gdp_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/all_B_to_GDP.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array_i
    close(unit = 20)

    do sim_ind = 1,N_sim_done
        tmp_array(sim_ind,:) = sim_res_all%SRA(sim_ind)%b_ex_to_GDP_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/all_B_ex_to_GDP.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array
    close(unit = 20)

    !Also - save all shock realizations (both indices and actual realizations) -
    !These will be useful when we compute various statistics
    do sim_ind = 1,N_sim_done
        tmp_array(sim_ind,:) = sim_res_all%SRA(sim_ind)%s_ind_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/all_s_ind.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array
    close(unit = 20)

    do sim_ind = 1,N_sim_done
        tmp_array_i(sim_ind,:,:) = sim_res_all%SRA(sim_ind)%g_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/all_g.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array_i
    close(unit = 20)

    !Also save consumption, labour, output, utility, and SWF (1-dimensional).
    do sim_ind = 1,N_sim_done
        tmp_array_i(sim_ind,:,:) = sim_res_all%SRA(sim_ind)%c_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/all_c.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array_i
    close(unit = 20)



    do sim_ind = 1,N_sim_done
        tmp_array_i(sim_ind,:,:) = sim_res_all%SRA(sim_ind)%l_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/all_l.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array_i
    close(unit = 20)

    do sim_ind = 1,N_sim_done
        tmp_array_i(sim_ind,:,:) = sim_res_all%SRA(sim_ind)%Y_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/all_Y.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array_i
    close(unit = 20)

    do sim_ind = 1,N_sim_done
        tmp_array_i(sim_ind,:,:) = sim_res_all%SRA(sim_ind)%U_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/all_u.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array_i
    close(unit = 20)

    do sim_ind = 1,N_sim_done
        tmp_array(sim_ind,:) = sim_res_all%SRA(sim_ind)%SWF_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/all_SWF.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array
    close(unit = 20)

    do sim_ind = 1,N_sim_done
        tmp_array(sim_ind,:) = sim_res_all%SRA(sim_ind)%term_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/all_term.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array
    close(unit = 20)

    !For debugging purposes, also save states a_prime, aex_prime, and rho_prime. We will then
    !plot these along with the grid boundaries (which will be read from the saved grids which are
    !one folder above the oher results). To some extent, this replicates what we already do
    !in Fortran (computing and saving a_prime_05, etc. But it is more informative because we see
    !more in the fan charts than just the extreme quantiles. Also it is worth preserving the stuff
    !in Fortran becuase it can be used quickly to detect cases where boundaries are too narrow,
    !without having to download the results from HPC, going to Windows to generate plots in Matlab, etc.).
    do sim_ind = 1,N_sim_done
        tmp_array_i(sim_ind,:,:) = sim_res_all%SRA(sim_ind)%a_prime_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/all_a_prime.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array_i
    close(unit = 20)

    do sim_ind = 1,N_sim_done
        tmp_array(sim_ind,:) = sim_res_all%SRA(sim_ind)%rho_prime_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/all_rho_prime.out'
    open(unit=20, file = file_path, status = 'replace')
    write(unit=20, fmt=*) tmp_array
    close(unit = 20)

    do sim_ind = 1,N_sim_done
        tmp_array(sim_ind,:) = sim_res_all%SRA(sim_ind)%aex_prime_sim
    end do
    file_path = 'results/'//trim(folder_name)//'/sim/all_aex_prime.out'
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

!Subroutine loss_func is used in subroutine LFFC (looking for feasible choice), when
!we are trying to find a point which minimizes the loss.
subroutine loss_func(n,xc,fc,iuser,ruser,inform)
        real (dp), intent (out) :: fc
        integer, intent (in) :: n
        real (dp), intent (inout) :: ruser(*)
        real (dp), intent (in) :: xc(n)
        integer, intent (inout) :: iuser(*)
        integer, intent(out) :: inform

        !pointer pars_pntr points to parameters type (it is more convenient to pass it using a common
        !block than create copies of everything we need and put them in a long common block). We
        !just need to be careful to not make any assignment statements because we would then modify
        !parameters elsewhere in the program!!!
        type(par), pointer :: pars_pntr

        !Values of state variables as saved in the grid (hence the suffix st).
        real(dp), dimension(I_par) :: a_st
        real(dp), dimension(I_par-1) :: rho_st
        real(dp) :: aex_st
        integer :: g_st_ind

        !Allocation and next-period states
        real(dp), dimension(M_par,I_par) :: c,l
        real(dp), dimension(M_par,I_par) :: a_pr
        real(dp), dimension(M_par,1) :: rho_pr, aex_pr

        logical :: constr_ok

        integer :: getavfail !IF this is greater than zero, then the choice xc implies a negative consumption or
        !labour supply of some agent in some state. In this case, it is hopeless to try to locally improve such allocation
        !set loss to max_loss and quit.

        !We need to pass pars to this function, and current states (the latter object varies with every call).
        !We need to be careful not to change anything in pars_pntr
        common /pars_pointer/ pars_pntr
        common /loss_func_bloc/ a_st,rho_st,aex_st,g_st_ind

        inform = 0

        !Compute all the remaining variables implied by the choice xc
        call get_av(pars_pntr,xc,a_st,rho_st,aex_st,g_st_ind,getavfail,c,l,a_pr,rho_pr,aex_pr,.false.)


        !If non-negativity of consumption or labour is violated, then current return is not defined. It is probably pointless
        !to try any local improvements. Just set inform = -1 which stops the maximization (and set loss to the worst possible
        !one)
        if(getavfail > 0) then
            fc = max_loss !maximum loss
            return
        end if

        !Compute the loss. I think that it is a good idea to use adjusted loss (so that we are looking for strictly feasible points
        !rather than simply points which violate the bounds as little as possible and could sit on the corner of
        !feasible set, very far from any optimum.
        call check_constr(pars_pntr,c,l,a_pr,rho_pr,aex_pr,constr_ok,fc,.false.,pars_pntr%LFFC_neg_loss)
end subroutine loss_func




!Subroutine fix_c is used to adjust a consumption vector such that it satisfies the minimum
!requirement that consumption is non-negative for all agents and all states. This is
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

!Subroutine get_lambda computes Lagrange multipliers lambda (one per each country)
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


!____TEMPORARY VARIABLES FOR DEBUGGING PURPOSES (remove later)____________________
!    real(dp) :: A_a,B_a,C_a
!    real(dp) :: Ec1recip,c1_imp_err,x
!    integer :: m
!_________________________________________________________________________________

    sing_fail = .false.

    !First pre-compute some useful expressions
    call util_c(c,u_c,pars%A_par,1.0_dp)
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

    !DEBUG - try regularizing before inversion (to see if the issue of near singularity is caused by
    !bad use of the subroutine or debug
    !TEST
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


!    write(*,*) (matmul(K_sys_copy,x_sys) - b_sys)
!    write(*,*) 'Stopping in subroutine get_lambda'
!    error stop

!    !TEST
!    A_a = pars%mass(1,1)*pars%alpha(1,1)
!    B_a = lambda(s_ind,1)*a_st(1)*pars%P(g_st_ind,s_ind)
!    C_a = lambda(s_ind,1)*(pars%B_par * l(s_ind,1)**(1.0_dp + pars%gamma_par) - 1.0_dp) !ok
!    !The following line is different in sign for agent 2 (if we ever want to do the test for changes in agent 2's consumption)
!    C_a = C_a - omega(s_ind)*( (1.0_dp/theta_lvl(1,1)) * pars%B_par * l(s_ind,1)**(pars%gamma_par) )
!    C_a = C_a - xi(s_ind)*pars%mass(1,1)
!    C_a = C_a - pars%beta**2.0_dp * R(s_ind,1) * a_prime(s_ind,1) * lambda_e(s_ind,1)
!
!    !Ec1recip
!    Ec1recip = 0.0_dp
!    do m = 1,M_par
!        Ec1recip = Ec1recip + (1.0_dp/c(m,1)) * pars%P(g_st_ind,m)
!    end do
!
!    x = c(s_ind,1)
!
!    !The error in equation (remember: x = c_pl_a(s_ind_a,1))
!    c1_imp_err = A_a/x + B_a * ((x*Ec1recip)**(-2.0_dp)) + C_a
!
!    write(*,*) 'c1_imp_err = ',c1_imp_err
    end do

end subroutine get_lambda

!(this function needs to be adapted for the new program).
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


end module mod_union

