Graetz Problem 2D Heat Transfer Solver
-------------------

Repo containing a simple c++ script for solving the Graetz Heat-Transfer Problem.

Graetz Problem: Temperatur distribution in fully-developed fluid flow of uniform inlet temperature flowing into a pipe of circular cross-section where the pipe walls are held at a constant temperature.

A good description of the Graetz oroblem can be found in the paper "The Graetz Problem" by R. Shankar Subramanian included here in the repository as 'Graetz_Problem_Subramanian.pdf' and available online [here](https://web2.clarkson.edu/projects/subramanian/ch490/notes/Graetz%20Problem.pdf).

Some other mathematical setup notes can be found in the file 'Graetz_Problem_Setup_Zachary_Cotman.pdf'

2 Versions of the program are available: one that solves the problem without including conduction, and one that includes conduction.

## Exectution of the makefile and scripts require:
* make
* g++
* gnu scientific libraries

## Files

* 'graetz_problem.cpp' - main file for convection-only solution
* 'graetz_problem_conduct.cpp' - main file for solution which includes conduction
* 'graetz_residual_functions.cpp' - functions for calculating the residuals (errors)
* 'graetz_xy_functions.cpp' - functions are used to build the right hand side(rhs), y, and initial guess, x, vectors used in the GSL sparse-matrix solver GMRES
* 'graetz_coeff_mat_functions.cpp' - functions are used to build the coefficient matricies used in the GSL sparse-matrix solver GMRES
* 'input_parameters.h' - boundary conditions and 2D mesh parameters, etc

* 'graetz_problem.dat' - output temperature values
* 'graetz_problem_conduct.dat' - output temperature values
* 'graetz_problem_residuals.dat' - output residuals
* 'graetz_problem.info' - some info from the simulation - useful primarily for debugging
* 'graetz_problem_conduct.info' - some info from the simulation - useful primarily for debugging

* 'make_graetz_problem' - makefile ( make -f make_graetz_problem )
* 'make_graetz_problem_conduct' - makefile ( make -f make_graetz_problem_conduct )

## Post Processing

Post-processing and plotting scripts originally made using Mathematica.
I may re-write them in Jupytern Notebooks to make this accessible.

See the file 'Graetz_Problem_Zachary_Cotman.pdf' for the original report containing some plots from the post-processing scripts. Note that I actually made a mistake in a simple calculation in the setup.
I may revisit, fix, and replcae this if I write some post-processing Jupyter Notebookes.

