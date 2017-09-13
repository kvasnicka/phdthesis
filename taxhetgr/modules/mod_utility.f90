module mod_utility
!This is the module containing utility function and other functions
!which depend on the functional form of utility function.
!
!
!When the functional form is changed, this module has to be replaced.
!This is computationally more efficient than branching during program execution
!and it is only slightly inconvenient if we want to compare specifications
!across utility functional forms. Also can use pre-processing and conditional
!compilation later.

!There are potential efficiency gains by approaching this differently
!but it's not a priority at this stage - get the program to work and then
!experiment with changes and see what they do to runtime... (It's not 
!clear how various compiler options woould perform).

use mod_types
use mod_parameters_hard

implicit none


contains

!These functions should work for arbitrary M (but can assume I=2)

!The functions which depend on utility but not on the functional form
!should be in different module - mod_taxhetgr or elsewhere.


!subroutine util(c,l,u) computes utility of all agents in all states.
!c is consumption, l is labour, u is utility. Element (m,i) of
!these matrices refers to state m agent i value.
subroutine util(c,l,u,A_par,B_par,Gamma_par)
    real(dp), dimension(:,:), intent(in) :: c,l !Usually dimension(M_par,I_par)
    real(dp), dimension(:,:), intent(out) :: u
    real(dp), intent(in) :: A_par,B_par,Gamma_par !parameters in the utility function

    !This could be written as an elemental function but doing it
    !as  a subroutine makes it easier to introduce heterogeneity.

    u = A_par*log(c) - B_par*(l**(1.0_dp+Gamma_par))/(1.0_dp+Gamma_par)
end subroutine util

!subroutine util_c computes marginal utility of consumption of all agents in all states.
!c is consumption, u is utility. Element (m,i) of
!these matrices refers to state m agent i value.
subroutine util_c(c,u_c,A_par)
    real(dp), dimension(:,:), intent(in) :: c !usually dimension will be (M_par,I_par)
    real(dp), dimension(:,:), intent(out) :: u_c
    real(dp), intent(in) :: A_par

    u_c = A_par*(c**(-1))
end subroutine util_c

!subroutine util_c computes marginal utility of labour of all agents in all states.
!l is consumption, u_l is MU of labour. Element (m,i) of
!these matrices refers to state m agent i value.
subroutine util_l(l,u_l,B_par,Gamma_par)
    real(dp), dimension(:,:), intent(in) :: l !usually dimension will be (M_par,I_par)
    real(dp), dimension(:,:), intent(out) :: u_l
    real(dp), intent(in) :: B_par,Gamma_par

    u_l = -B_par * (l**Gamma_par)
end subroutine util_l

!The formulas used in the next couple of subroutines are given in the paper (appendix).

!subroutine lab returns labour supply of all agents in all states
!given consumption of all agents in all states.
!Input
!   - c: dim(M_par,I_par) - consumption of all agents in all states
!   - theta: dim(I_par) - productivity of all agents (levels) - one-dimensional becuase not stochastic
!       !and not used in matrix multiplication anywhere.
!   - g: government expenditure (in levels) for all shock realizations - as column vector
!   - mass: dim(1,I_par) - mass of agents (it's a row vector so we can use it in
!       matrix multiplication elsewhere.
!   - Gamma_par: parameter in utility function.
!Output
!   - l: dim(M_par,I_par) - labour supply of all agents in all states
subroutine lab(c,l,theta,mass,g,Gamma_par)
    real(dp), dimension(M_par,I_par), intent(in) :: c
    real(dp), dimension(M_par,I_par), intent(out) :: l
    real(dp), dimension(I_par), intent(in) ::theta
    real(dp), dimension(M_par,1), intent(in) :: g
    real(dp), dimension(1,I_par), intent(in) :: mass
    real(dp), intent(in) :: Gamma_par

    real(dp), dimension(M_par,1) :: tot_exp !consumption + government expenditure

    tot_exp(:,1) = mass(1,1)*c(:,1) + mass(1,2)*c(:,2) + g(:,1)
    !agent 1
    l(:,1) = tot_exp(:,1)/&
    (mass(1,1)*theta(1) + mass(1,2)*theta(2)*((theta(2)*c(:,1)/(theta(1)*c(:,2)))**(1.0_dp/Gamma_par)))
    !agent 2
    l(:,2) = tot_exp(:,1)/&
    (mass(1,1)*theta(1)*((theta(1)*c(:,2)/(theta(2)*c(:,1)))**(1.0_dp/Gamma_par)) +  mass(1,2)*theta(2))
end subroutine lab

!Subroutine lab_fp computes labour in the first period (it may be better to unify the subroutines
!lab and lab_fp but lab will need generalizing for I>2 anyway so can do it then).
subroutine lab_fp(c,l,theta,mass,g,Gamma_par)
    real(dp), dimension(1,I_par), intent(in) :: c
    real(dp), dimension(1,I_par), intent(out) :: l
    real(dp), dimension(1,I_par), intent(in) ::theta
    real(dp), intent(in) :: g !now a scalar
    real(dp), dimension(1,I_par), intent(in) :: mass
    real(dp), intent(in) :: Gamma_par

    real(dp) :: tot_exp !consumption + government expenditure

    tot_exp = mass(1,1)*c(1,1) + mass(1,2)*c(1,2) + g


    !agent 1
    l(1,1) = tot_exp/&
    (mass(1,1)*theta(1,1) + mass(1,2)*theta(1,2)*((theta(1,2)*c(1,1)/(theta(1,1)*c(1,2)))**(1.0_dp/Gamma_par)))
    !agent 2
    l(1,2) = tot_exp/&
    (mass(1,1)*theta(1,1)*((theta(1,1)*c(1,2)/(theta(1,2)*c(1,1)))**(1.0_dp/Gamma_par)) +  mass(1,2)*theta(1,2))
end subroutine lab_fp

!subroutine ass computes next-period marginal-utility-adjusted asset holdings
!for all agents in all states.
!Input:
!   - a_st: current-period marginal-utility-adjusted asset holdings (state variable) in levels
!(so we need to multiply by productivity first and can't just use the value read from the grid)
!   -   c,l: consumption and labour supply of all agents in all states
!   - A_par,B_par,Gamma_par, beta: parameters in the utility function and the discount factor
!   - P_onerow: row s_ind of transition matrix of the Markov process governing shocks
!       stacked on top of each other M times. This is useful for computing conditional expectation
!       (more efficient than do cycle) - expectation is computed conditional on the
!       current shock being with index s_ind. This is pre-compured in pars%P_onerow(s_ind,:,:).
!Output:
!   - a_pr : next-period marginal-utility-adjusted asset holdings (pr stands for 'prime', std. notation)
subroutine ass(a_lvl,c,l,A_par,B_par,Gamma_par,beta,P_onerow,a_pr)
    real(dp), dimension(I_par), intent(in) :: a_lvl
    real(dp), dimension(M_par,I_par), intent(in) :: c,l
    real(dp), intent(in) :: A_par,B_par,Gamma_par,beta
    real(dp), dimension(M_par,M_par), intent(in) :: P_onerow

    real(dp), dimension(M_par,I_par) :: a_pr

    !asset holdings for agent 1
    a_pr(:,1) = (a_lvl(1)/(c(:,1)* matmul(P_onerow,c(:,1)**(-1))) + B_par*(l(:,1)**(1.0_dp+gamma_par)) - A_par)/beta
    !agent 2
    a_pr(:,2) = (a_lvl(2)/(c(:,2)* matmul(P_onerow,c(:,2)**(-1))) + B_par*(l(:,2)**(1.0_dp+gamma_par)) - A_par)/beta
end subroutine ass

!Subroutine b_prime_fp computes the next-period asset holdings of every agent
!in the first period. This is not particularly efficient but not important
!because it's not called many times.
!
!It can be generalized so that it calls subroutines which compute
!marginal utility of c and l and then it wouldn't need to be modified
!when the utility function is modified.
subroutine b_prime_fp(b_init,c_vec,l_vec,A_par,B_par,Gamma_par,b_prime)
    real(dp), dimension(1,I_par), intent(in) :: b_init,c_vec,l_vec
    real(dp), intent(in) :: A_par,B_par,Gamma_par
    real(dp), dimension(1,I_par), intent(out) :: b_prime

    b_prime = b_init + c_vec*((B_par/A_par)*(l_vec**(1.0_dp+Gamma_par)) -1.0_dp)
end subroutine b_prime_fp

!function c2M_int (int stands for interior) returns agent 2's consumption in state M (the last one)
!given consumption of all agents in all states except for agent 2's consumption in state M,
!and a few other inputs. If some constraint on asset holding binds,
! iterative procedure is used instead of this function.
!
!This function depends on the functional form of the utility function.
!Inputs:
!   - c_inc: cons of all agents in all states except for agent 2's state M consumption (inc stands for incomplete)
!   - rho_st: ratio of last-period marginal utilities (agent2/agent1), state variable.
!   - P: transition matrix of the Markov chain shock process
!   - s_ind_st: the index of last-period shock realization (state variable)
function c2M_int(c_inc,rho_st,P,s_ind_st)
    real(dp) :: c2M_int

    real(dp), dimension(M_par,I_par), intent(in) :: c_inc
    real(dp), intent(in) :: rho_st
    real(dp), dimension(M_par,M_par), intent(in) :: P
    integer, intent(in) :: s_ind_st

    real(dp), dimension(1,1) :: tmpmat !for temporary storage of result


    tmpmat = rho_st*P(s_ind_st,M_par) / (matmul(P(s_ind_st:s_ind_st,:),c_inc(:,1:1)**(-1)) &
    - rho_st*matmul(P(s_ind_st:s_ind_st,1:M_par-1),c_inc(1:M_par-1,2:2)**(-1)) )

    c2M_int = tmpmat(1,1)

!(I > 2 generalization needed)

end function c2M_int


!Function c1_imp_err computes error as function of implied consumption c1_imp. This is error in the non-linear equation
!in c1(st), and if this error is zero, the equation is solved. If we want to compute this error for agent 2, then we
!have to change this slightly. c1_imp and lambda_t1_imp are scalars, corresponding to one state only.
!
!Warning: this function already presumes the logarithmic utility function. We will need to generalize this if
!we are interested in accuracy tests for other functions.
function c1_imp_err2(x,iuser,ruser)
    real(dp) :: c1_imp_err2
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

    c1_imp_err2 = 0.0_dp !initialize value

    !in c_pl, change the element (corresponding to agent 1's consumption in state s_ind) to the implied consumption
    !for which we are trying to compute the error
    c_pl_a(s_ind_a,1) = x

    !Ec1recip
    Ec1recip = 0.0_dp
    do m = 1,M_par
        Ec1recip = Ec1recip + (1.0_dp/c_pl_a(m,1)) * P_a(s_st_ind_a,m)
    end do

    !The error in equation (remember: x = c_pl_a(s_ind_a,1))
    c1_imp_err2 = A_a/x + B_a * ((x*Ec1recip)**(-2.0_dp)) + C_a

end function c1_imp_err2


end module mod_utility
