module mod_parameters_hard
!This module contains all the parameters of interest that are hardcoded
!and not contained in parameter file.

use mod_types

implicit none


!The following two variables should conceptually be in the file with 
!other parameters but it would create a lot of unnecessary difficulties
!and they will be changed only rarely.
integer, parameter :: I_par = 2 !Number of agents

integer, parameter :: M_par = 2 !Number of realizations of shock space.

integer, parameter :: ncvar_par = M_par*I_par - 1 !Number of choice variables

integer :: use_threads_mkl = 1 !Number of threads to be used by Intel Math Kernel Library.
!This should be 1 because we are using all of the physical cores to distribute workload
!on the value function. Also, MKL apparently uses OpenMP which doesn't work with CAF applications
!on Intel compiler.

integer :: spec_img = 1 !index of special image which will handle some synchronization, I/O, etc.
!Could perhaps achieve marginal improvement of performance if some other than 1 is used
!(for example the last image will usually have fewer operations in V iteration to do
!than the other images).


integer, parameter :: N_How = 3 !Number of thresholds in Howard acceleration algorithm implementation.
!(hardcoded for simplicity)

end module mod_parameters_hard
