module mod_utility
!This modle contains the utility function and other functions
!which depend on the functional form of utility function.
!
!
!When the functional form is changed, this module has to be replaced.
!This is computationally more efficient than branching during program execution
!and it is only slightly inconvenient if we want to compare specifications
!across utility functional forms.
!
!If we want to change a functional form of utility functions, we should
!only need to replace this module. However, in practice, there are a few
!other places in the code in which the functional form is assumed,
!mostly in module mod_union.f90. If we are interested in experimenting
!with different functional forms, this needs to be fixed in the future!
!
!The function in this module work for arbitrary number of shock realizations, but some of them
!they assume that there are two countries (I_par = 2). We need to generalize this
!if we want to solve the model for more countries in the future.

use mod_types
use mod_parameters_hard

implicit none

contains

!Subroutine a_prime computes next-period marginal-adjusted debt of both countries.
!(for all shock realizations). As input, we need the current asset holdings (a_st), consumption and labour
!supply of all countries in all states (c,l), pre-computed expected marginal utilities of consumption in both
!countries (Eu_c), and some parameters of utility function.
!As input, we need the current as
subroutine a_prime(a_st,c,l,Eu_c,mass,theta,G,B_par,sigma_par,Gamma_par,beta,a_pr)
    real(dp), dimension(I_par), intent(in) :: a_st
    real(dp), dimension(M_par,I_par), intent(in) :: c,l
    real(dp), dimension(1,I_par), intent(in) :: Eu_c !this is explicitly a row vector so we do not have to use reshape
    !after computing this using matrix multiplication elsewhere.
    real(dp), dimension(1,I_par), intent(in) :: mass,theta
    real(dp), dimension(M_par,I_par), intent(in) :: G !government expenditure in level (column i contains expenditures in all states in both countries)

    real(dp), intent(in) :: B_par,sigma_par,Gamma_par,beta

    real(dp), dimension(M_par,I_par) :: a_pr

    !asset holdings for country 1
    a_pr(:,1) = (c(:,1)**(-sigma_par)/beta) * (a_st(1)/Eu_c(1,1) + mass(1,1) * &
    ( G(:,1) - theta(1,1)*l(:,1) + B_par * c(:,1)**(sigma_par) * l(:,1)**(1.0_dp + gamma_par) ) )

    !asset holdings for country 2
    a_pr(:,2) = (c(:,2)**(-sigma_par)/beta) * (a_st(2)/Eu_c(1,2) + mass(1,2) * &
    ( G(:,2)  - theta(1,2)*l(:,2) + B_par * c(:,2)**(sigma_par) * l(:,2)**(1.0_dp + gamma_par) ) )

end subroutine a_prime

!Subroutine aex_prime computes next-period marginal-adjusted external debt of country 1.
!(for all shock realizations). Uses a subset of the inputs used by a_prime
subroutine aex_prime(aex_st,c,l,Eu_c,mass,theta,G,sigma_par,beta,aex_pr)
    real(dp), intent(in) :: aex_st
    real(dp), dimension(M_par,I_par), intent(in) :: c,l
    real(dp), dimension(1,I_par), intent(in) :: Eu_c !this is explicitly a row vector so we do not have to use reshape
    !after computing this using matrix multiplication elsewhere.
    real(dp), dimension(1,I_par), intent(in) :: mass,theta
    real(dp), dimension(M_par,I_par), intent(in) :: G !government expenditure in level (column i contains expenditures in all states in both countries)

    real(dp), intent(in) :: sigma_par,beta

    real(dp), dimension(M_par,1) :: aex_pr !next-period external debt of country one for all shock realizations

    aex_pr(:,1) = (c(:,1)**(-sigma_par)/beta) * (aex_st/Eu_c(1,1) + mass(1,1) * (G(:,1) + c(:,1) - theta(1,1)*l(:,1) ) )
end subroutine aex_prime

!subroutine b_ex_prime_fp does the same but for external debt
subroutine b_ex_prime_fp(b_ex_init,c,l,G_lvl,mass,theta_0,b_ex_prime)
    real(dp), dimension(1,I_par), intent(in) :: c,l,G_lvl,mass,theta_0
    real(dp), intent(in) :: b_ex_init

    real(dp), dimension(1,1), intent(out) :: b_ex_prime !the dimension just for consistency with code elsewhere

    b_ex_prime(1,1) = b_ex_init + mass(1,1)*(G_lvl(1,1) + c(1,1) - theta_0(1,1)*l(1,1))

end subroutine b_ex_prime_fp

!Subroutine b_prime_fp computes the next-period debt of every country
!in the first period
subroutine b_prime_fp(b_init,c,l,G_lvl,mass,theta_0,B_par,sigma_par,gamma_par,b_prime)
    real(dp), dimension(1,I_par), intent(in) :: b_init,c,l,G_lvl,mass,theta_0
    real(dp), intent(in) :: B_par,sigma_par,gamma_par

    real(dp), dimension(1,I_par), intent(out) :: b_prime

    b_prime = b_init + mass*G_lvl - mass * (1.0_dp - B_par * (c**sigma_par) * ((l**gamma_par)/theta_0) )*theta_0 * l

end subroutine b_prime_fp

!Subroutine c_IM computes consumption in country I in shock realization M,
!given the other elements of the consumption matrix, and the state variable rho.
!This function depends on the functional form of the utility function.
!Inputs:
!   - c: cons of all agents in all states except for agent 2's state M consumption
!   - rho_st: ratio of last-period marginal utilities (agent2/agent1), state variable.
!   - P: transition matrix of the Markov chain shock process
!   - s_ind_st: the index of last-period shock realization (state variable)
subroutine c_IM(c,rho_st,sigma,P,s_ind_st)

    real(dp), dimension(M_par,I_par), intent(inout) :: c !On entry - the incomplete consumption matrix (missing element M,I)
    !on exit it will be completed.
    real(dp), intent(in) :: rho_st
    real(dp), intent(in) :: sigma !parameter in utility function
    real(dp), dimension(M_par,M_par), intent(in) :: P
    integer, intent(in) :: s_ind_st

    if(sigma == 1.0_dp) then
    c(M_par:M_par,I_par:I_par) = rho_st*P(s_ind_st,M_par) / (matmul(P(s_ind_st:s_ind_st,:),c(:,1:1)**(-sigma)) &
    - rho_st*matmul(P(s_ind_st:s_ind_st,1:M_par-1),c(1:M_par-1,2:2)**(-sigma)) )
    else
    c(M_par:M_par,I_par:I_par) = (rho_st*P(s_ind_st,M_par) / (matmul(P(s_ind_st:s_ind_st,:),c(:,1:1)**(-sigma)) &
    - rho_st*matmul(P(s_ind_st:s_ind_st,1:M_par-1),c(1:M_par-1,2:2)**(-sigma)) ))**(1.0_dp/sigma)
    end if

!(I > 2 generalization needed)
end subroutine c_IM

!a debug version of the subroutine. It displays some output.
subroutine c_IM_debug(c,rho_st,sigma,P,s_ind_st)

    real(dp), dimension(M_par,I_par), intent(inout) :: c !On entry - the incomplete consumption matrix (missing element M,I)
    !on exit it will be completed.
    real(dp), intent(in) :: rho_st
    real(dp), intent(in) :: sigma !parameter in utility function
    real(dp), dimension(M_par,M_par), intent(in) :: P
    integer, intent(in) :: s_ind_st

    if(sigma == 1.0_dp) then
    c(M_par:M_par,I_par:I_par) = rho_st*P(s_ind_st,M_par) / (matmul(P(s_ind_st:s_ind_st,:),c(:,1:1)**(-sigma)) &
    - rho_st*matmul(P(s_ind_st:s_ind_st,1:M_par-1),c(1:M_par-1,2:2)**(-sigma)) )
    else
    c(M_par:M_par,I_par:I_par) = (rho_st*P(s_ind_st,M_par) / (matmul(P(s_ind_st:s_ind_st,:),c(:,1:1)**(-sigma)) &
    - rho_st*matmul(P(s_ind_st:s_ind_st,1:M_par-1),c(1:M_par-1,2:2)**(-sigma)) ))**(1.0_dp/sigma)
    end if

write(*,*) '~~~~~~~~~Debug in subroutine c_IM~~~~~~~~~~~~'
    write(*,*) 'rho_st = ',rho_st
    write(*,*) 'c = ',c
    write(*,*) 'P = ',P
    write(*,*) 'sigma = ',sigma
    write(*,*) 's_ind_st = ',s_ind_st
    !write(*,*)

write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

!(I > 2 generalization needed)
end subroutine c_IM_debug

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

    !(This subroutine is redundant and should be removed eventually).
    if(this_image() == 1) then
        write(*,*) 'WARNING: Subroutine lab is used even though it corresponds to the old program.'
    end if

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


!subroutine util(c,l,u) computes utility of all agents in all states.
!c is consumption, l is labour, u is utility. Element (m,i) of
!these matrices refers to state m agent i value.
subroutine util(c,l,u,A_par,B_par,sigma,Gamma_par)
    real(dp), dimension(:,:), intent(in) :: c,l !Usually dimension(M_par,I_par)
    real(dp), dimension(:,:), intent(out) :: u
    real(dp), intent(in) :: A_par,B_par,sigma,Gamma_par !parameters in the utility function

    if(sigma == 1.0_dp) then !(logarithmic utility)
        u = A_par*log(c) - B_par*(l**(1.0_dp+Gamma_par))/(1.0_dp+Gamma_par)
    else
        u = A_par*(c**(1.0_dp-sigma))/(1.0_dp-sigma) - B_par*(l**(1.0_dp+Gamma_par))/(1.0_dp+Gamma_par)
    end if

end subroutine util

!subroutine util_c computes marginal utility of consumption of all agents in all states.
!c is consumption, u is utility. Element (m,i) of
!these matrices refers to state m agent i value.
subroutine util_c(c,u_c,A,sigma)
    real(dp), dimension(:,:), intent(in) :: c !usually dimension will be (M_par,I_par)
    real(dp), dimension(:,:), intent(out) :: u_c
    real(dp), intent(in) :: A,sigma !parameters in the utility function

    u_c = A*(c**(-sigma))
end subroutine util_c

!subroutine util_l computes marginal utility of labour of all agents in all states.
!l is consumption, u_l is MU of labour. Element (m,i) of
!these matrices refers to state m agent i value.
subroutine util_l(l,u_l,B_par,Gamma_par)
    real(dp), dimension(:,:), intent(in) :: l !usually dimension will be (M_par,I_par)
    real(dp), dimension(:,:), intent(out) :: u_l
    real(dp), intent(in) :: B_par,Gamma_par

    u_l = -B_par * (l**Gamma_par)
end subroutine util_l

end module mod_utility
