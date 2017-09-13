This repository contains the programs for solving the models in my PhD thesis

Folder taxhetgr contains the program for solving the model of Chapter 2 (Heterogeneous labour productivity
trends and optimal fiscal policy)

Folder taxhetgrtr contains the program for solving the model of Chapter 3 (Heterogeneous labour productivity
trends and optimal fiscal policy with transfers)

Folder union contains the program for solving the model of Chapter 4 (Optimal fiscal policy in a currency union)

Each of the three folders contains a folder 'parameters', in which there is a sample parameter file (baseline.txt),
which contains all relevant parameters. The name of the parameter file can be passed to the program via a command line
argument. For example, if the executable file is 'union', running the command 'union baseline' loads the parameters
from file baseline.txt in folder parameters. If the program is run without argument, the parameter file
parameters/baseline.txt is used.


A high-level overview of the algorithms is in the PhD thesis. The codes themselves contain more detailed comments
where appropriate.

The programs are written in Fortran (2008 standard) and require a compiler which supports Coarray Fortran.

The codes use NAG Fortran numerical library (Mark 24), and Intel MKL libraries.
