#Makefile for the program taxhetgr (coarray using shared memory)
#Some of the options vary depending on whether we compile on local computer or on HPC
#(need to comment/uncomment appropriate options)

#Variables which control output
#CAF is the option for coarray fortran. Use either coarray=shared or coarray=distributed
#For use on multiple clusters on HPC, use -coarray=distributed
CAF = -coarray=shared
#CAF = -coarray=distributed

#Libraries should run on one core only because the paralelization happens at higher
#level
MKL = -mkl=sequential

#_______Definition of NAG_related_variables_________________
#Add $(NAG_INC) to the end of every ifort command which compiles something which uses NAG libraries (MKL version)

#Path to NAG library on local system (default used in make all)
#Local version
NAG_DIR = /opt/NAG/fll6i24dcl
#HPC version (also need to load module for that)
#NAG_DIR = /usr/local/Cluster-Apps/nag/fll6i23dcl

NAG_INC_mkl =  -I $(NAG_DIR)/nag_interface_blocks $(NAG_DIR)/lib/libnag_mkl.a -Wl,--start-group $(NAG_DIR)/mkl_intel64/libmkl_intel_lp64.a $(NAG_DIR)/mkl_intel64/libmkl_intel_thread.a $(NAG_DIR)/mkl_intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread
NAG_INC_nag = -I $(NAG_DIR)/nag_interface_blocks $(NAG_DIR)/lib/libnag_nag.a

#Here, choose either NAG_INC_mkl or NAG_INC_NAG (the latter being slower but potentially more reliable)
NAG_INC = $(NAG_INC_mkl)
#NAG_INC = $(NAG_INC_nag)

#______________End of definition of NAG related variables_____________

#Optimization options massively increase compilation time so for local development/debugging
#use OPT_flag = and only use optimization such as OPT_flag = -O3 -xhost for when we're submitting
#computations on HPC or long-running computations locally
#No optimization (above standard)
OPT_flag =
OPT_flag = -O3 -xhost
#Options such as check all and check bounds are useful but make the runtime really slow.

#debug_option
#debug_option = -g -debug all 
#debug_option = -g -check bounds
debug_option = -g

#Heap arrays option (to avoid segmentation faults caused by large arrays)
#Do not use - leads to errors when program is run on multiple nodes and
#shared arrays are large!!! Instead use sys_call_result =
#system('ulimit -Ss unlimited')
HA =
#before had HA = -heap-arrays 10

#To do (at some point, low priority)
#define variables containing multiple files (e.g., mod_types.o, mod_par.o, ... to reduce length of the lines)

all: mod_types.o mod_IO.o mod_par.o mod_parameters_hard.o mod_tools.o mod_work_dist.o mod_union.o mod_utility.o bspline_sub_module.o union.o
	ifort $(OPT_flag) $(CAF) $(MKL) $(HA) $(debug_option) -o union mod_types.o mod_IO.o mod_par.o mod_parameters_hard.o mod_tools.o mod_work_dist.o mod_union.o  mod_utility.o bspline_sub_module.o union.o  $(NAG_INC)
		
mod_types.o: modules/mod_types.f90
	ifort $(OPT_flag) $(HA) -c $(debug_option) modules/mod_types.f90 $(NAG_INC)
	
mod_IO.o: mod_types.o modules/mod_IO.f90
	ifort $(OPT_flag) $(HA) -c $(debug_option) modules/mod_IO.f90 
	
mod_par.o: mod_types.o mod_parameters_hard.o mod_IO.o modules/mod_par.f90
	ifort $(OPT_flag) $(HA) -c $(debug_option) modules/mod_par.f90
	
mod_parameters_hard.o: mod_types.o modules/mod_parameters_hard.f90
	ifort $(OPT_flag) $(HA) -c $(debug_option) modules/mod_parameters_hard.f90
	
union.o: mod_types.o mod_IO.o mod_par.o mod_parameters_hard.o mod_tools.o mod_work_dist.o union.f90
	ifort $(OPT_flag) $(CAF) $(MKL) $(HA) -c -C $(debug_option) union.f90 $(NAG_INC)
	
mod_tools.o: mod_types.o mod_IO.o modules/mod_tools.f90
	ifort $(OPT_flag) $(CAF) $(HA) -c $(debug_option) modules/mod_tools.f90
		
mod_work_dist.o: mod_types.o mod_union.o modules/mod_work_dist.f90
	ifort $(OPT_flag) $(CAF) $(HA) -c $(debug_option) modules/mod_work_dist.f90
	
mod_union.o: mod_utility.o mod_types.o mod_par.o mod_parameters_hard.o mod_tools.o bspline_sub_module.o modules/mod_union.f90
	ifort $(OPT_flag) $(CAF) $(HA) -c $(debug_option) modules/mod_union.f90 $(NAG_INC) 
	
mod_utility.o: mod_types.o mod_parameters_hard.o modules/mod_utility.f90
	ifort $(OPT_flag) $(CAF) $(HA) -c $(debug_option)  modules/mod_utility.f90
	
bspline_sub_module.o: modules/bspline_sub_module.f90
	ifort $(OPT_flag) $(CAF) $(HA) -c $(debug_option)  modules/bspline_sub_module.f90
		
#coarray_debugging_mod.o: modules/coarray_debugging_mod.f90
#	ifort $(OPT_flag) $(CAF) $(HA) -c $(debug_option)  modules/coarray_debugging_mod.f90
	

#clean deletes the executable and all *.o and *.mod files in the folder
clean:
	rm -f union *.o *.mod
	
	
#End of makefile
