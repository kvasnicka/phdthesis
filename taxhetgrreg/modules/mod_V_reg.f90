module mod_V_reg
!(file mod_V_reg.f90 automatically generated by subroutine gen_mod_V_reg)
 
!Module mod_V_reg contains hard-coded constants and functional expressions which determine
!how the value function is approximated. This module is the only thing that needs to be
!changed when we change the set of explanatory variables (used to approximate value function).
 
use mod_types
use mod_parameters_hard
 
implicit none
 
!Definitions of variables
 
!V_exp_num is the number of explanatory variables in the regressions (including the constant term
!if it is present and excluding the (M_par - 1) shock realization dummies which are always present).
!This has to correspond to the actual number of variables defined in the functions/
!subroutines below. V_exp_num_dyn is the number of variables in the dynamic case.
integer, parameter :: V_exp_num = 216
integer, parameter :: V_exp_num_dyn = 5
!V_reg_basis is the GAF_basis in parameter file which was used to generate this file.
!We save this here because the value in parameter file may have changed since generating
!this module, and particularly since this module was compiled. Therefore, using the value
!in the parameter file could lead to errors. Some parts of the program depend on what basis
!was used. For example, if we are using Cheb. pol., then states are normalized into 0,1 interval.
 
 integer, dimension(2), parameter :: V_reg_basis = [3,2]
 
!(GEN: in the later project with 4 value function, there will be 4 separate constants V1_exp_num,...,V4_exp_num).
 
contains
 
!function V_und_exp transforms the underlying states (those that appear in the recursive formulation)
!into their functions which serve as explanatory variables in the regressions approximating the value
!function. V_und_exp_dyn is the version for dynamic case.
!(GEN: in the later project, there will be 4 subroutines: V1_und_exp,...,V4_und_exp)
!To help avoid mistakes, the inputs of this subroutine suggest which state they represent (consistent
!with notation used in the recursive formulation).
function V_und_exp(a1,a2,rho)
   real(dp), dimension(V_exp_num) :: V_und_exp
   real(dp), intent(in) :: a1
   real(dp), intent(in) :: a2
   real(dp), intent(in) :: rho
   !states explanatory in regression (functions of underlying states):



 
   !Tensor product of Cheb. pol. (pars%GAF_basis = 3)
   !The function assumes that the inputs were already projected onto [-1,1] interval.
   V_und_exp(1) = 1.0_dp !constant term
   V_und_exp(2) = a1   !a1
   V_und_exp(3) = 2*V_und_exp(2)*a1-V_und_exp(1)   !a1**2
   V_und_exp(4) = 2*V_und_exp(3)*a1-V_und_exp(2)   !a1**3
   V_und_exp(5) = 2*V_und_exp(4)*a1-V_und_exp(3)   !a1**4
   V_und_exp(6) = 2*V_und_exp(5)*a1-V_und_exp(4)   !a1**5
   V_und_exp(7) = a2   !a2
   V_und_exp(8) = a1*a2   !a1*a2
   V_und_exp(9) = V_und_exp(3)*a2   !a1**2*a2
   V_und_exp(10) = V_und_exp(4)*a2   !a1**3*a2
   V_und_exp(11) = V_und_exp(5)*a2   !a1**4*a2
   V_und_exp(12) = V_und_exp(6)*a2   !a1**5*a2
   V_und_exp(13) = 2*V_und_exp(7)*a2-V_und_exp(1)   !a2**2
   V_und_exp(14) = a1*V_und_exp(13)   !a1*a2**2
   V_und_exp(15) = V_und_exp(3)*V_und_exp(13)   !a1**2*a2**2
   V_und_exp(16) = V_und_exp(4)*V_und_exp(13)   !a1**3*a2**2
   V_und_exp(17) = V_und_exp(5)*V_und_exp(13)   !a1**4*a2**2
   V_und_exp(18) = V_und_exp(6)*V_und_exp(13)   !a1**5*a2**2
   V_und_exp(19) = 2*V_und_exp(13)*a2-V_und_exp(7)   !a2**3
   V_und_exp(20) = a1*V_und_exp(19)   !a1*a2**3
   V_und_exp(21) = V_und_exp(3)*V_und_exp(19)   !a1**2*a2**3
   V_und_exp(22) = V_und_exp(4)*V_und_exp(19)   !a1**3*a2**3
   V_und_exp(23) = V_und_exp(5)*V_und_exp(19)   !a1**4*a2**3
   V_und_exp(24) = V_und_exp(6)*V_und_exp(19)   !a1**5*a2**3
   V_und_exp(25) = 2*V_und_exp(19)*a2-V_und_exp(13)   !a2**4
   V_und_exp(26) = a1*V_und_exp(25)   !a1*a2**4
   V_und_exp(27) = V_und_exp(3)*V_und_exp(25)   !a1**2*a2**4
   V_und_exp(28) = V_und_exp(4)*V_und_exp(25)   !a1**3*a2**4
   V_und_exp(29) = V_und_exp(5)*V_und_exp(25)   !a1**4*a2**4
   V_und_exp(30) = V_und_exp(6)*V_und_exp(25)   !a1**5*a2**4
   V_und_exp(31) = 2*V_und_exp(25)*a2-V_und_exp(19)   !a2**5
   V_und_exp(32) = a1*V_und_exp(31)   !a1*a2**5
   V_und_exp(33) = V_und_exp(3)*V_und_exp(31)   !a1**2*a2**5
   V_und_exp(34) = V_und_exp(4)*V_und_exp(31)   !a1**3*a2**5
   V_und_exp(35) = V_und_exp(5)*V_und_exp(31)   !a1**4*a2**5
   V_und_exp(36) = V_und_exp(6)*V_und_exp(31)   !a1**5*a2**5
   V_und_exp(37) = rho   !rho
   V_und_exp(38) = a1*rho   !a1*rho
   V_und_exp(39) = V_und_exp(3)*rho   !a1**2*rho
   V_und_exp(40) = V_und_exp(4)*rho   !a1**3*rho
   V_und_exp(41) = V_und_exp(5)*rho   !a1**4*rho
   V_und_exp(42) = V_und_exp(6)*rho   !a1**5*rho
   V_und_exp(43) = a2*rho   !a2*rho
   V_und_exp(44) = a1*a2*rho   !a1*a2*rho
   V_und_exp(45) = V_und_exp(3)*a2*rho   !a1**2*a2*rho
   V_und_exp(46) = V_und_exp(4)*a2*rho   !a1**3*a2*rho
   V_und_exp(47) = V_und_exp(5)*a2*rho   !a1**4*a2*rho
   V_und_exp(48) = V_und_exp(6)*a2*rho   !a1**5*a2*rho
   V_und_exp(49) = V_und_exp(13)*rho   !a2**2*rho
   V_und_exp(50) = a1*V_und_exp(13)*rho   !a1*a2**2*rho
   V_und_exp(51) = V_und_exp(3)*V_und_exp(13)*rho   !a1**2*a2**2*rho
   V_und_exp(52) = V_und_exp(4)*V_und_exp(13)*rho   !a1**3*a2**2*rho
   V_und_exp(53) = V_und_exp(5)*V_und_exp(13)*rho   !a1**4*a2**2*rho
   V_und_exp(54) = V_und_exp(6)*V_und_exp(13)*rho   !a1**5*a2**2*rho
   V_und_exp(55) = V_und_exp(19)*rho   !a2**3*rho
   V_und_exp(56) = a1*V_und_exp(19)*rho   !a1*a2**3*rho
   V_und_exp(57) = V_und_exp(3)*V_und_exp(19)*rho   !a1**2*a2**3*rho
   V_und_exp(58) = V_und_exp(4)*V_und_exp(19)*rho   !a1**3*a2**3*rho
   V_und_exp(59) = V_und_exp(5)*V_und_exp(19)*rho   !a1**4*a2**3*rho
   V_und_exp(60) = V_und_exp(6)*V_und_exp(19)*rho   !a1**5*a2**3*rho
   V_und_exp(61) = V_und_exp(25)*rho   !a2**4*rho
   V_und_exp(62) = a1*V_und_exp(25)*rho   !a1*a2**4*rho
   V_und_exp(63) = V_und_exp(3)*V_und_exp(25)*rho   !a1**2*a2**4*rho
   V_und_exp(64) = V_und_exp(4)*V_und_exp(25)*rho   !a1**3*a2**4*rho
   V_und_exp(65) = V_und_exp(5)*V_und_exp(25)*rho   !a1**4*a2**4*rho
   V_und_exp(66) = V_und_exp(6)*V_und_exp(25)*rho   !a1**5*a2**4*rho
   V_und_exp(67) = V_und_exp(31)*rho   !a2**5*rho
   V_und_exp(68) = a1*V_und_exp(31)*rho   !a1*a2**5*rho
   V_und_exp(69) = V_und_exp(3)*V_und_exp(31)*rho   !a1**2*a2**5*rho
   V_und_exp(70) = V_und_exp(4)*V_und_exp(31)*rho   !a1**3*a2**5*rho
   V_und_exp(71) = V_und_exp(5)*V_und_exp(31)*rho   !a1**4*a2**5*rho
   V_und_exp(72) = V_und_exp(6)*V_und_exp(31)*rho   !a1**5*a2**5*rho
   V_und_exp(73) = 2*V_und_exp(37)*rho-V_und_exp(1)   !rho**2
   V_und_exp(74) = a1*V_und_exp(73)   !a1*rho**2
   V_und_exp(75) = V_und_exp(3)*V_und_exp(73)   !a1**2*rho**2
   V_und_exp(76) = V_und_exp(4)*V_und_exp(73)   !a1**3*rho**2
   V_und_exp(77) = V_und_exp(5)*V_und_exp(73)   !a1**4*rho**2
   V_und_exp(78) = V_und_exp(6)*V_und_exp(73)   !a1**5*rho**2
   V_und_exp(79) = a2*V_und_exp(73)   !a2*rho**2
   V_und_exp(80) = a1*a2*V_und_exp(73)   !a1*a2*rho**2
   V_und_exp(81) = V_und_exp(3)*a2*V_und_exp(73)   !a1**2*a2*rho**2
   V_und_exp(82) = V_und_exp(4)*a2*V_und_exp(73)   !a1**3*a2*rho**2
   V_und_exp(83) = V_und_exp(5)*a2*V_und_exp(73)   !a1**4*a2*rho**2
   V_und_exp(84) = V_und_exp(6)*a2*V_und_exp(73)   !a1**5*a2*rho**2
   V_und_exp(85) = V_und_exp(13)*V_und_exp(73)   !a2**2*rho**2
   V_und_exp(86) = a1*V_und_exp(13)*V_und_exp(73)   !a1*a2**2*rho**2
   V_und_exp(87) = V_und_exp(3)*V_und_exp(13)*V_und_exp(73)   !a1**2*a2**2*rho**2
   V_und_exp(88) = V_und_exp(4)*V_und_exp(13)*V_und_exp(73)   !a1**3*a2**2*rho**2
   V_und_exp(89) = V_und_exp(5)*V_und_exp(13)*V_und_exp(73)   !a1**4*a2**2*rho**2
   V_und_exp(90) = V_und_exp(6)*V_und_exp(13)*V_und_exp(73)   !a1**5*a2**2*rho**2
   V_und_exp(91) = V_und_exp(19)*V_und_exp(73)   !a2**3*rho**2
   V_und_exp(92) = a1*V_und_exp(19)*V_und_exp(73)   !a1*a2**3*rho**2
   V_und_exp(93) = V_und_exp(3)*V_und_exp(19)*V_und_exp(73)   !a1**2*a2**3*rho**2
   V_und_exp(94) = V_und_exp(4)*V_und_exp(19)*V_und_exp(73)   !a1**3*a2**3*rho**2
   V_und_exp(95) = V_und_exp(5)*V_und_exp(19)*V_und_exp(73)   !a1**4*a2**3*rho**2
   V_und_exp(96) = V_und_exp(6)*V_und_exp(19)*V_und_exp(73)   !a1**5*a2**3*rho**2
   V_und_exp(97) = V_und_exp(25)*V_und_exp(73)   !a2**4*rho**2
   V_und_exp(98) = a1*V_und_exp(25)*V_und_exp(73)   !a1*a2**4*rho**2
   V_und_exp(99) = V_und_exp(3)*V_und_exp(25)*V_und_exp(73)   !a1**2*a2**4*rho**2
   V_und_exp(100) = V_und_exp(4)*V_und_exp(25)*V_und_exp(73)   !a1**3*a2**4*rho**2
   V_und_exp(101) = V_und_exp(5)*V_und_exp(25)*V_und_exp(73)   !a1**4*a2**4*rho**2
   V_und_exp(102) = V_und_exp(6)*V_und_exp(25)*V_und_exp(73)   !a1**5*a2**4*rho**2
   V_und_exp(103) = V_und_exp(31)*V_und_exp(73)   !a2**5*rho**2
   V_und_exp(104) = a1*V_und_exp(31)*V_und_exp(73)   !a1*a2**5*rho**2
   V_und_exp(105) = V_und_exp(3)*V_und_exp(31)*V_und_exp(73)   !a1**2*a2**5*rho**2
   V_und_exp(106) = V_und_exp(4)*V_und_exp(31)*V_und_exp(73)   !a1**3*a2**5*rho**2
   V_und_exp(107) = V_und_exp(5)*V_und_exp(31)*V_und_exp(73)   !a1**4*a2**5*rho**2
   V_und_exp(108) = V_und_exp(6)*V_und_exp(31)*V_und_exp(73)   !a1**5*a2**5*rho**2
   V_und_exp(109) = 2*V_und_exp(73)*rho-V_und_exp(37)   !rho**3
   V_und_exp(110) = a1*V_und_exp(109)   !a1*rho**3
   V_und_exp(111) = V_und_exp(3)*V_und_exp(109)   !a1**2*rho**3
   V_und_exp(112) = V_und_exp(4)*V_und_exp(109)   !a1**3*rho**3
   V_und_exp(113) = V_und_exp(5)*V_und_exp(109)   !a1**4*rho**3
   V_und_exp(114) = V_und_exp(6)*V_und_exp(109)   !a1**5*rho**3
   V_und_exp(115) = a2*V_und_exp(109)   !a2*rho**3
   V_und_exp(116) = a1*a2*V_und_exp(109)   !a1*a2*rho**3
   V_und_exp(117) = V_und_exp(3)*a2*V_und_exp(109)   !a1**2*a2*rho**3
   V_und_exp(118) = V_und_exp(4)*a2*V_und_exp(109)   !a1**3*a2*rho**3
   V_und_exp(119) = V_und_exp(5)*a2*V_und_exp(109)   !a1**4*a2*rho**3
   V_und_exp(120) = V_und_exp(6)*a2*V_und_exp(109)   !a1**5*a2*rho**3
   V_und_exp(121) = V_und_exp(13)*V_und_exp(109)   !a2**2*rho**3
   V_und_exp(122) = a1*V_und_exp(13)*V_und_exp(109)   !a1*a2**2*rho**3
   V_und_exp(123) = V_und_exp(3)*V_und_exp(13)*V_und_exp(109)   !a1**2*a2**2*rho**3
   V_und_exp(124) = V_und_exp(4)*V_und_exp(13)*V_und_exp(109)   !a1**3*a2**2*rho**3
   V_und_exp(125) = V_und_exp(5)*V_und_exp(13)*V_und_exp(109)   !a1**4*a2**2*rho**3
   V_und_exp(126) = V_und_exp(6)*V_und_exp(13)*V_und_exp(109)   !a1**5*a2**2*rho**3
   V_und_exp(127) = V_und_exp(19)*V_und_exp(109)   !a2**3*rho**3
   V_und_exp(128) = a1*V_und_exp(19)*V_und_exp(109)   !a1*a2**3*rho**3
   V_und_exp(129) = V_und_exp(3)*V_und_exp(19)*V_und_exp(109)   !a1**2*a2**3*rho**3
   V_und_exp(130) = V_und_exp(4)*V_und_exp(19)*V_und_exp(109)   !a1**3*a2**3*rho**3
   V_und_exp(131) = V_und_exp(5)*V_und_exp(19)*V_und_exp(109)   !a1**4*a2**3*rho**3
   V_und_exp(132) = V_und_exp(6)*V_und_exp(19)*V_und_exp(109)   !a1**5*a2**3*rho**3
   V_und_exp(133) = V_und_exp(25)*V_und_exp(109)   !a2**4*rho**3
   V_und_exp(134) = a1*V_und_exp(25)*V_und_exp(109)   !a1*a2**4*rho**3
   V_und_exp(135) = V_und_exp(3)*V_und_exp(25)*V_und_exp(109)   !a1**2*a2**4*rho**3
   V_und_exp(136) = V_und_exp(4)*V_und_exp(25)*V_und_exp(109)   !a1**3*a2**4*rho**3
   V_und_exp(137) = V_und_exp(5)*V_und_exp(25)*V_und_exp(109)   !a1**4*a2**4*rho**3
   V_und_exp(138) = V_und_exp(6)*V_und_exp(25)*V_und_exp(109)   !a1**5*a2**4*rho**3
   V_und_exp(139) = V_und_exp(31)*V_und_exp(109)   !a2**5*rho**3
   V_und_exp(140) = a1*V_und_exp(31)*V_und_exp(109)   !a1*a2**5*rho**3
   V_und_exp(141) = V_und_exp(3)*V_und_exp(31)*V_und_exp(109)   !a1**2*a2**5*rho**3
   V_und_exp(142) = V_und_exp(4)*V_und_exp(31)*V_und_exp(109)   !a1**3*a2**5*rho**3
   V_und_exp(143) = V_und_exp(5)*V_und_exp(31)*V_und_exp(109)   !a1**4*a2**5*rho**3
   V_und_exp(144) = V_und_exp(6)*V_und_exp(31)*V_und_exp(109)   !a1**5*a2**5*rho**3
   V_und_exp(145) = 2*V_und_exp(109)*rho-V_und_exp(73)   !rho**4
   V_und_exp(146) = a1*V_und_exp(145)   !a1*rho**4
   V_und_exp(147) = V_und_exp(3)*V_und_exp(145)   !a1**2*rho**4
   V_und_exp(148) = V_und_exp(4)*V_und_exp(145)   !a1**3*rho**4
   V_und_exp(149) = V_und_exp(5)*V_und_exp(145)   !a1**4*rho**4
   V_und_exp(150) = V_und_exp(6)*V_und_exp(145)   !a1**5*rho**4
   V_und_exp(151) = a2*V_und_exp(145)   !a2*rho**4
   V_und_exp(152) = a1*a2*V_und_exp(145)   !a1*a2*rho**4
   V_und_exp(153) = V_und_exp(3)*a2*V_und_exp(145)   !a1**2*a2*rho**4
   V_und_exp(154) = V_und_exp(4)*a2*V_und_exp(145)   !a1**3*a2*rho**4
   V_und_exp(155) = V_und_exp(5)*a2*V_und_exp(145)   !a1**4*a2*rho**4
   V_und_exp(156) = V_und_exp(6)*a2*V_und_exp(145)   !a1**5*a2*rho**4
   V_und_exp(157) = V_und_exp(13)*V_und_exp(145)   !a2**2*rho**4
   V_und_exp(158) = a1*V_und_exp(13)*V_und_exp(145)   !a1*a2**2*rho**4
   V_und_exp(159) = V_und_exp(3)*V_und_exp(13)*V_und_exp(145)   !a1**2*a2**2*rho**4
   V_und_exp(160) = V_und_exp(4)*V_und_exp(13)*V_und_exp(145)   !a1**3*a2**2*rho**4
   V_und_exp(161) = V_und_exp(5)*V_und_exp(13)*V_und_exp(145)   !a1**4*a2**2*rho**4
   V_und_exp(162) = V_und_exp(6)*V_und_exp(13)*V_und_exp(145)   !a1**5*a2**2*rho**4
   V_und_exp(163) = V_und_exp(19)*V_und_exp(145)   !a2**3*rho**4
   V_und_exp(164) = a1*V_und_exp(19)*V_und_exp(145)   !a1*a2**3*rho**4
   V_und_exp(165) = V_und_exp(3)*V_und_exp(19)*V_und_exp(145)   !a1**2*a2**3*rho**4
   V_und_exp(166) = V_und_exp(4)*V_und_exp(19)*V_und_exp(145)   !a1**3*a2**3*rho**4
   V_und_exp(167) = V_und_exp(5)*V_und_exp(19)*V_und_exp(145)   !a1**4*a2**3*rho**4
   V_und_exp(168) = V_und_exp(6)*V_und_exp(19)*V_und_exp(145)   !a1**5*a2**3*rho**4
   V_und_exp(169) = V_und_exp(25)*V_und_exp(145)   !a2**4*rho**4
   V_und_exp(170) = a1*V_und_exp(25)*V_und_exp(145)   !a1*a2**4*rho**4
   V_und_exp(171) = V_und_exp(3)*V_und_exp(25)*V_und_exp(145)   !a1**2*a2**4*rho**4
   V_und_exp(172) = V_und_exp(4)*V_und_exp(25)*V_und_exp(145)   !a1**3*a2**4*rho**4
   V_und_exp(173) = V_und_exp(5)*V_und_exp(25)*V_und_exp(145)   !a1**4*a2**4*rho**4
   V_und_exp(174) = V_und_exp(6)*V_und_exp(25)*V_und_exp(145)   !a1**5*a2**4*rho**4
   V_und_exp(175) = V_und_exp(31)*V_und_exp(145)   !a2**5*rho**4
   V_und_exp(176) = a1*V_und_exp(31)*V_und_exp(145)   !a1*a2**5*rho**4
   V_und_exp(177) = V_und_exp(3)*V_und_exp(31)*V_und_exp(145)   !a1**2*a2**5*rho**4
   V_und_exp(178) = V_und_exp(4)*V_und_exp(31)*V_und_exp(145)   !a1**3*a2**5*rho**4
   V_und_exp(179) = V_und_exp(5)*V_und_exp(31)*V_und_exp(145)   !a1**4*a2**5*rho**4
   V_und_exp(180) = V_und_exp(6)*V_und_exp(31)*V_und_exp(145)   !a1**5*a2**5*rho**4
   V_und_exp(181) = 2*V_und_exp(145)*rho-V_und_exp(109)   !rho**5
   V_und_exp(182) = a1*V_und_exp(181)   !a1*rho**5
   V_und_exp(183) = V_und_exp(3)*V_und_exp(181)   !a1**2*rho**5
   V_und_exp(184) = V_und_exp(4)*V_und_exp(181)   !a1**3*rho**5
   V_und_exp(185) = V_und_exp(5)*V_und_exp(181)   !a1**4*rho**5
   V_und_exp(186) = V_und_exp(6)*V_und_exp(181)   !a1**5*rho**5
   V_und_exp(187) = a2*V_und_exp(181)   !a2*rho**5
   V_und_exp(188) = a1*a2*V_und_exp(181)   !a1*a2*rho**5
   V_und_exp(189) = V_und_exp(3)*a2*V_und_exp(181)   !a1**2*a2*rho**5
   V_und_exp(190) = V_und_exp(4)*a2*V_und_exp(181)   !a1**3*a2*rho**5
   V_und_exp(191) = V_und_exp(5)*a2*V_und_exp(181)   !a1**4*a2*rho**5
   V_und_exp(192) = V_und_exp(6)*a2*V_und_exp(181)   !a1**5*a2*rho**5
   V_und_exp(193) = V_und_exp(13)*V_und_exp(181)   !a2**2*rho**5
   V_und_exp(194) = a1*V_und_exp(13)*V_und_exp(181)   !a1*a2**2*rho**5
   V_und_exp(195) = V_und_exp(3)*V_und_exp(13)*V_und_exp(181)   !a1**2*a2**2*rho**5
   V_und_exp(196) = V_und_exp(4)*V_und_exp(13)*V_und_exp(181)   !a1**3*a2**2*rho**5
   V_und_exp(197) = V_und_exp(5)*V_und_exp(13)*V_und_exp(181)   !a1**4*a2**2*rho**5
   V_und_exp(198) = V_und_exp(6)*V_und_exp(13)*V_und_exp(181)   !a1**5*a2**2*rho**5
   V_und_exp(199) = V_und_exp(19)*V_und_exp(181)   !a2**3*rho**5
   V_und_exp(200) = a1*V_und_exp(19)*V_und_exp(181)   !a1*a2**3*rho**5
   V_und_exp(201) = V_und_exp(3)*V_und_exp(19)*V_und_exp(181)   !a1**2*a2**3*rho**5
   V_und_exp(202) = V_und_exp(4)*V_und_exp(19)*V_und_exp(181)   !a1**3*a2**3*rho**5
   V_und_exp(203) = V_und_exp(5)*V_und_exp(19)*V_und_exp(181)   !a1**4*a2**3*rho**5
   V_und_exp(204) = V_und_exp(6)*V_und_exp(19)*V_und_exp(181)   !a1**5*a2**3*rho**5
   V_und_exp(205) = V_und_exp(25)*V_und_exp(181)   !a2**4*rho**5
   V_und_exp(206) = a1*V_und_exp(25)*V_und_exp(181)   !a1*a2**4*rho**5
   V_und_exp(207) = V_und_exp(3)*V_und_exp(25)*V_und_exp(181)   !a1**2*a2**4*rho**5
   V_und_exp(208) = V_und_exp(4)*V_und_exp(25)*V_und_exp(181)   !a1**3*a2**4*rho**5
   V_und_exp(209) = V_und_exp(5)*V_und_exp(25)*V_und_exp(181)   !a1**4*a2**4*rho**5
   V_und_exp(210) = V_und_exp(6)*V_und_exp(25)*V_und_exp(181)   !a1**5*a2**4*rho**5
   V_und_exp(211) = V_und_exp(31)*V_und_exp(181)   !a2**5*rho**5
   V_und_exp(212) = a1*V_und_exp(31)*V_und_exp(181)   !a1*a2**5*rho**5
   V_und_exp(213) = V_und_exp(3)*V_und_exp(31)*V_und_exp(181)   !a1**2*a2**5*rho**5
   V_und_exp(214) = V_und_exp(4)*V_und_exp(31)*V_und_exp(181)   !a1**3*a2**5*rho**5
   V_und_exp(215) = V_und_exp(5)*V_und_exp(31)*V_und_exp(181)   !a1**4*a2**5*rho**5
   V_und_exp(216) = V_und_exp(6)*V_und_exp(31)*V_und_exp(181)   !a1**5*a2**5*rho**5
 
end function V_und_exp
 
function V_und_exp_dyn(a1,a2,rho,t)
   real(dp), dimension(V_exp_num_dyn) :: V_und_exp_dyn
   real(dp), intent(in) :: a1
   real(dp), intent(in) :: a2
   real(dp), intent(in) :: rho
   real(dp), intent(in) :: t
   !states explanatory in regression (functions of underlying states):
 
   !Complete polynomials (pars%GAF_basis = 2)
   !K (sum of degreees) is the max. indiv. degree in pars%GAF_deg
   V_und_exp_dyn(1) = 1.0_dp !constant term
   V_und_exp_dyn(2) = a1   !a1
   V_und_exp_dyn(3) = a2   !a2
   V_und_exp_dyn(4) = rho   !rho
   V_und_exp_dyn(5) = t   !t
 
end function V_und_exp_dyn
 
 
!The following functions are used to scale features.
!This is extremely important for speed of convergence in iterative algorithms
!like gradient descent. If use_FS is false (which is set in parameter file,
!then these functions do not change anything.
 
function V_und_exp_FS(st_exp,use_FS)
   real(dp), dimension(V_exp_num) :: V_und_exp_FS
   real(dp), dimension(V_exp_num), intent(in) :: st_exp !vector of exp vars before scaling
   logical :: use_FS !If false, no scaling is done.
 
       V_und_exp_FS = st_exp
       return
 
end function V_und_exp_FS
 
function V_und_exp_dyn_FS(st_exp,use_FS)
   real(dp), dimension(V_exp_num_dyn) :: V_und_exp_dyn_FS
   real(dp), dimension(V_exp_num_dyn), intent(in) :: st_exp !vector of exp vars before scaling
   logical :: use_FS !If false, no scaling is done.
 
   if(.not. use_FS) then
       V_und_exp_dyn_FS = st_exp
       return
   end if

   V_und_exp_dyn_FS(1) = (st_exp(1))
   V_und_exp_dyn_FS(2) = (st_exp(2) - 0.000_dp)/6.00000000000000_dp
   V_und_exp_dyn_FS(3) = (st_exp(3) - 0.000_dp)/6.00000000000000_dp
   V_und_exp_dyn_FS(4) = (st_exp(4) - 0.240_dp)/0.120000000000000_dp
   V_und_exp_dyn_FS(5) = (st_exp(5) - 0.000_dp)/1.000000000000000E-013_dp
 
end function V_und_exp_dyn_FS
 
 
!Note: no interpolation over shock realizations (discrete Markov chain). Dummy variables introduced for
!each possible realization (so we essentially treat these as separate functions).
 
!Function V_eval evaluates the approximated value function. Inputs are coefficients (V_coeff) in the
!regression approximation as column vector, and a value of states (obtained as output of subroutine V_und_exp).
!These inputs need to be of the same dimension (no check done for performance considerations - this function is
!very deep in the program).
 
function V_eval(V_coeff,st_exp)
   real(dp) :: V_eval !experiment
   real(dp), dimension(:), intent(in) :: st_exp !generated by function V_und_exp
   real(dp), dimension(:,:), intent(in) :: V_coeff !coefficients stored in column vector
 
   V_eval = dot_product(st_exp,reshape(V_coeff,shape(st_exp)))
end function V_eval
end module mod_V_reg
