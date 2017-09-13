module mod_parameters_hard
!This module contains all the parameters of interest that are hardcoded
!and not contained in parameter file.

use mod_types

implicit none

!The following two variables should conceptually be in the file with 
!other parameters but it would create a lot of unnecessary difficulties
!and they will be changed only rarely.
integer, parameter :: I_par = 2 !Number of agents

integer, parameter :: M_par = 2 !Number of realizations of shock. If shocks are i.i.d. (which
!is set in the parameter file), then M_par is not used in grid construction.

integer, parameter :: ncvar_par = M_par*(I_par+1) - 1 !Number of choice variables


integer :: use_threads_mkl = 1 !Number of threads to be used by Intel Math Kernel Library.
!This should be 1 because we are using all of the physical cores to distribute workload
!on the value function. Also, MKL apparently uses OpenMP which doesn't work with CAF applications
!on Intel compiler.

integer, parameter :: spec_img = 1 !index of special image which will handle some synchronization, I/O, etc.
!This should not be change from the default value of 1 at this stage (other images than 1 sometimes
!have issues with i/o, etc.)

integer, parameter :: N_How = 3 !Number of thresholds in Howard acceleration algorithm implementation.
!(hardcoded for simplicity)

real(dp), parameter :: max_loss = 1000000000000.0_dp !Maximum value of loss function (if this is exceeded, there
!is no point knowing the exact value anyway, or attempting local improvements.

real(dp), parameter, dimension(M_par,1) :: ones_col = 1.0_dp!column vector of ones

real(dp), parameter, dimension(I_par,1) :: ones_col_I = 1.0_dp!column vector of ones - dimension I (number of countries)


end module mod_parameters_hard
