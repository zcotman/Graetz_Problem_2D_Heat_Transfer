
//***********************************************************************************
//  FILE: graetz_problem.cpp
// 
//  Heat Transfer
//  
//  Determining the Temperature Distribution
//      With Conduction
//		Pipe of Circular Cross-Section
//		Fully Developed Flow
//		BC:	Tin(x=0,r) known
//			Twall(x,r = r0) constant 
//
//  Programmer:  Zach Cotman  cotman.6@osu.edu
//
//  Revision history:
//      Nov 20, 2018:	Original Version
//		Nov 21, 2018:	Using some funcitons for readability
//		Nov 22, 2018:	Added header file for inputs 
//
//  Notes:
//   Adapted from the GSL Sparse Matrix Tutorial at
//		https://www.gnu.org/software/gsl/doc/html/splinalg.html
//	 A big thank-you to 
//	 	Dr. Richard J Furnstahl and
//	 	Dr. Ralf A Bundschuh
//    for teaching me everything I know about c++
//
//***********************************************************************************

// include files
#include <iostream>
#include <iomanip> 
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <string>
using namespace std; // meh

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>

#include "graetz_problem_conduct.h"
// Input Parameters Which Adjust the Mesh and Solution Space 
#include "input_parameters.h"

//___________________________________________________________________________________
// MAIN

int
main ()
{
	//*******************************************************************************
	// SOME SETUP
	//*******************************************************************************
		// Output File Stream for Simulation Info
		ofstream out_h;
		out_h.open("graetz_problem_conduct.info");

		// beta
		const double beta = alpha;
		
		out_h 	<< "beta   = " << beta   << "\n"
				<< "radius = " << radius << "\n"
				<< "length = " << length << "\n"
				<< "dx     = " << dx     << "\n"
				<< "dr     = " << dr     << "\n"
				<< "N      = " << N      << "\n"
				<< "M      = " << M      << "\n";
				
		// some counters
		int i = 0, j = 0;
		
		//***************************************************************************
		// r, u, initial guess for temperature distribution

		gsl_vector   *r = gsl_vector_alloc(M);		// r position "vector"
		gsl_vector   *u = gsl_vector_alloc(M);		// fluid velocity "vector"
		gsl_matrix   *T = gsl_matrix_alloc(M, N);

		// Set r-position vector
		for(i = 0; i < M; ++i){gsl_vector_set(r, i, i*dr);}
	
		// Set velocity distribution vector
		// u = 1-r^2 - could read in from input file
		gsl_vector_set_all(u, -1.);
		gsl_vector_mul(u, r);
		gsl_vector_mul(u, r);
		gsl_vector_add_constant(u, 1.);
	
		// Input File Stream for Initial Temperature Distribution
		ifstream in;
		
		// Open "Marching" Data to use as initial guess for iterating program
		in.open("graetz_problem.dat");
		
		// get T data from file
		double data_point = 0.;
		for(i = 0; i < N; ++i){
			for(j = 0; j < M; ++j){
				in >> data_point;
				gsl_matrix_set(T, j, i, data_point);
			}
		}
		in.close(); // Close input file stream for initial temperature distribution
	
	// Allocate space for Matricies and Vectors using GSL Library
	
		// Coefficient Matricies
		gsl_spmatrix *Coeff_Interior = gsl_spmatrix_alloc(M, M);// "triplet format"
		gsl_spmatrix *Coeff_Interior_Compressed;				// "compressed format"
	
		gsl_spmatrix *Coeff_Final = gsl_spmatrix_alloc(M, M); 	// "triplet format"
		gsl_spmatrix *Coeff_Final_Compressed;					// "compressed format"
		
		const double  A = 1./dr/dr;				// A coefficient        - see notes
		const double  B = beta/dx/dx;			// B coefficient        - see notes
		gsl_vector   *C = gsl_vector_alloc(M);	// C coefficirnt vector - see notes
		gsl_vector   *D = gsl_vector_alloc(M);	// D coefficirnt vector - see notes
	
		gsl_vector   *y = gsl_vector_alloc(M);  // "rhs vector"
		gsl_vector   *x = gsl_vector_alloc(M);	// "solution vector"
	
	// Assign Elements to the Vectors C and D
	
		// C = 1/(r*2*dr)
		for(i = 0; i < M; ++i){gsl_vector_set(C, i, 1./2./i/dr/dr);}
	
		// D = u_j/(2dx)
		gsl_vector_memcpy(D, u);
		gsl_vector_scale(D, 1./dx);
	
		out_h << "B = " << B << "   A = " << A << "\n";
	
		out_h << "C = " << "\n";
		for(i = 0; i < M; ++i){
			out_h << gsl_vector_get(C,i) << " ";
		}
		out_h << "\n";
	
		out_h << "D = " << "\n";
		for(i = 0; i < M; ++i){
			out_h << gsl_vector_get(D,i) << " ";
		}
		out_h << "\n";
	
	//*******************************************************************************
	// COEFFICIENT MATRICIES
	//*******************************************************************************
	
	// Interior Columns
		set_coeff_interior(M, A, B, C, D, Coeff_Interior);
		// "convert to compressed column format"
		Coeff_Interior_Compressed = gsl_spmatrix_ccs(Coeff_Interior);
		// free some workspace
		gsl_spmatrix_free(Coeff_Interior);		
	
	// Final Column
		set_coeff_final(M, A, B, C, D, Coeff_Final);
		// "convert to compressed column format"
		Coeff_Final_Compressed = gsl_spmatrix_ccs(Coeff_Final);	
		// free workspace
		gsl_spmatrix_free(Coeff_Final);
		
	// My Solver Setup
		//int my_max_iter = 10000;
		int my_min_iter = 10;
		//The total number of iterations allowed is 100
		gsl_vector *Residuals = gsl_vector_alloc(my_max_iter);
		gsl_vector_set_all(Residuals, 0.);
		gsl_vector *Slopes = gsl_vector_alloc(my_max_iter);
		gsl_vector_set_all(Slopes, 0.);
		double res_tot = 0.; // the residual norm for a step
		int my_counter = 0; // number of iterations of my solver
		double slope = 1;	// change in residual
		double my_tol = 0.001;	// end criterion
		double delta = 0.;

	// GMRES: initiate and allocate memory for solver
		const double tol = 1.0e-4;	// "solution relative tolerance"
		const size_t max_iter = 10; // "maximum iterations"
		const gsl_splinalg_itersolve_type *S = gsl_splinalg_itersolve_gmres;
		gsl_splinalg_itersolve *work = gsl_splinalg_itersolve_alloc(S, M, 10);
		size_t iter = 0;
		int status; // will return the status of the iterative solver employed by GSL
	
	//*******************************************************************************
	// Run Solver
	//*******************************************************************************
	
	do{	// an iteration
		
	
	// update temperatures column by column
		
		// Step Through Interior Columns
		for(i = 1; i < N-1; ++i){
		
			// define current rhs vector "y"
			set_rhs_interior(M, i, B, D, T, y);
			
			// initial guess for solution vector
			set_initial_guess(M, i, x, T);
		
			// Solve Column using GSL GMRES
				// iterate and return status
				iter = 0;
				do {	
					status = gsl_splinalg_itersolve_iterate
					(Coeff_Interior_Compressed, y, tol, x, work);
				}
				while(status ==  GSL_CONTINUE && ++iter < max_iter);
		
			// Check to see if the GSL routine is iterating a lot
			// if(iter == max_iter)
			//     { cout << "i = " << i << " Max iterations reached \n"; }
		
			// Assign to T Matrix
			update_T(M, i, x, T);
		}
		
		// Final COlumn
		
			i = N-1;
			// set rhs vector "y"
			set_rhs_final(M, i, B, D, T, y);
			
			// initial guess for solution vector x
			set_initial_guess(M, i, x, T);
		
			// Solve Column using GSL GMRES
				// iterate and return status
				iter = 0;
				do {
					status = gsl_splinalg_itersolve_iterate
					(Coeff_Final_Compressed, y, tol, x, work);
				}
				while(status ==  GSL_CONTINUE && ++iter < max_iter);
		
			// Check to see if the GSL routine is iterating a lot
			// if(iter == max_iter)
			//     { cout << "i = " << i << " Max iterations reached \n"; }
		
			// Assign to T Matrix
			update_T(M, i, x, T);
		
		// calculate residuals
			
			// interior nodes
			res_tot += res_interior(N, M, A, B, C, D, T); 	// see external functions
			// "bottom" nodes
			res_tot += res_bottom(N, M, A, B, C, D, T);		// see external functions
			// last line of nodes
			res_tot += res_lastline(N, M, A, B, C, D, T);	// see external functions
			// bottom right node
			res_tot += res_corner(N, M, A, B, C, D, T);		// see external functions
			
			// total residual
			res_tot = sqrt(res_tot)/(N-1)/(M-1);
		
			// save in residuals vector
			gsl_vector_set(Residuals, my_counter, res_tot);
			
			// calculate slope of residuals change
			if(my_counter > 1){
				slope = (gsl_vector_get(Residuals, my_counter-2) 
				      - gsl_vector_get(Residuals, my_counter))/2.;
				gsl_vector_set(Slopes, my_counter - 2, slope);
				delta = slope/gsl_vector_get(Slopes, 0);
			}
			
			++my_counter;
			
			if(my_counter % 100 == 0){cout << my_counter << "\n";}
			
	} // End of "do" loop
	// Stop itterating if the residuals are not changing by much 
	//	(relative to the maximum residual)
	while( my_counter < my_max_iter 
		/*|| (delta > my_tol &&  my_counter < my_max_iter)*/ );
		
	if(my_counter == my_max_iter){ cout << " MY Max iterations reached \n"; }

	//*******************************************************************************
	// Output
	//*******************************************************************************
	
	// Residuals
		ofstream out_r;									// open file stream
		out_r.open("graetz_problem_residuals.dat");		// residuals file
		for(i = 0; i < my_counter-1; ++i){
			out_r << gsl_vector_get(Residuals, i) << " ";
		}
		out_r.close();									// close file stream
	
	// Temperature Distribution
		ofstream out_d;								// declare output file stream
		out_d.open("graetz_problem_conduct.dat");	
		for(i = 0; i < N; ++i){
			for(j = 0; j < M; ++j){
				out_d << gsl_matrix_get(T, j, i) << " ";
			}
			out_d << "\n";
		}
		out_d.close();								// close file stream
	
	
	// Clean Up
	
	// free workspace allocated by gsl
	gsl_splinalg_itersolve_free(work);
	
	gsl_spmatrix_free(Coeff_Interior_Compressed);
	gsl_spmatrix_free(Coeff_Final_Compressed);
	
	gsl_vector_free(r);
	gsl_vector_free(u);
	gsl_vector_free(C);
	gsl_vector_free(D);
	gsl_vector_free(x);		
	gsl_vector_free(y);
	
	gsl_matrix_free(T);
	
	//Close Header File Stream
	out_h.close();
	
	return 0;
}//end of main
