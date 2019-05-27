
//***********************************************************************************
//  FILE: graetz_problem.cpp
// 
//  Heat Transfer
//  
//  Determining the Temperature Distribution
//      Convection Only
//		Pipe of Circular Cross-Section
//		Fully Developed Flow
//		BC:	Tin(x=0,r) known
//			Twall(x,r = r0) constant 
//
//  Programmer:  Zach Cotman  cotman.6@osu.edu
//
//  Revision history:
//      Nov 11, 2018:	Original Version
//		Nov 19, 2018:	Cleaned up some of the code
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

// Input Parameters Which Adjust the Mesh and Solution Space
#include "input_parameters.h"

//___________________________________________________________________________________
// MAIN

int
main ()
{
	//*******************************************************************************
	// Output File Stream
	ofstream out_d;						// declare output file stream
	out_d.open("graetz_problem.dat");
	ofstream out_h;
	out_h.open("graetz_problem.info");  // information about the simulation
	
	ofstream out;
	out.open("dx_dr.dat");
	out << dx << " " << dr;
	out.close();
	
	const double beta = 0.;
	
	out_h 	<< "Beta   = " << beta   << "\n"
			<< "radius = " << radius << "\n"
			<< "length = " << length << "\n"
			<< "dx     = " << dx     << "\n"
			<< "dr     = " << dr     << "\n"
			<< "N      = " << N      << "\n"
			<< "M      = " << M      << "\n";
			
	gsl_vector *Tin = gsl_vector_alloc(M); // allocate memory for input temperature
	gsl_vector_set_all(Tin, 1.);  // set input temperature
	gsl_vector_set(Tin, M-1, 0.); // wall temperature
	
	gsl_vector   *r = gsl_vector_alloc(M);	// r position "vector"
	gsl_vector   *u = gsl_vector_alloc(M);	// fluid velocity "vector"
	
	// r
	for(int i = 0; i < M; ++i){gsl_vector_set(r, i, i*dr);}
	
	// u = 1-r^2
	gsl_vector_set_all(u, -1.);
	gsl_vector_mul(u, r);
	gsl_vector_mul(u, r);
	gsl_vector_add_constant(u, 1.);
	
	// Ouptut for checking
	out_h << "Tin = " << "\n";
	for(int i = 0; i < M; ++i){
		out_h << gsl_vector_get(Tin,i) << " ";
	}
	out_h << "\n";
	
	out_h << "r = " << "\n";
	for(int i = 0; i < M; ++i){
		out_h << gsl_vector_get(r,i) << " ";
	}
	out_h << "\n";
	
	out_h << "u = " << "\n";
	for(int i = 0; i < M; ++i){
		out_h << gsl_vector_get(u,i) << " ";
	}
	out_h << "\n";

	

	//*******************************************************************************
	// Allocate space for Matricies and Vectors using GSL Library
	gsl_spmatrix *Coeff_Mat = gsl_spmatrix_alloc(M, M); // "triplet format"
	gsl_spmatrix *Coeff_Mat_Compressed;					// "compressed format"
	
	const double  A = 1./dr/dr;					// A coefficitne        - see notes
	const double  B = beta/dx/dx;				// B coefficitne        - see notes
	gsl_vector   *C = gsl_vector_alloc(M);		// C coefficitne vector - see notes
	gsl_vector   *D = gsl_vector_alloc(M);		// D coefficitne vector - see notes
	
	gsl_vector   *y = gsl_vector_alloc(M);      // "rhs vector"
	gsl_vector   *x = gsl_vector_alloc(M);		// "solution vector"
	
	//*******************************************************************************
	// Assign Elements to the Vectors C and D
	//	and the coefficient matrix Coeff_Mat
	
	// C = 1/(r*2*dr)
	for(int i = 0; i < M; ++i){gsl_vector_set(C, i, 1./2./i/dr/dr);}
	
	// D = u_j/(2dx)
	gsl_vector_memcpy(D, u);
	gsl_vector_scale(D, 1./2./dx);
	
	// output for checking
	out_h << "A = " << A << "\n";
	out_h << "B = " << B << "\n";

	out_h << "C = " << "\n";
	for(int i = 0; i < M; ++i){
		out_h << gsl_vector_get(C,i) << " ";
	}
	out_h << "\n";
	
	out_h << "D = " << "\n";
	for(int i = 0; i < M; ++i){
		out_h << gsl_vector_get(D,i) << " ";
	}
	out_h << "\n";

	//*******************************************************************************
	//construc Coeff_Mat
	double Cj = 0, Dj = 0;
	
	//out_h << "(Cj-A)/(B+Dj)/2.    (Dj+A+B)/(B+Dj)     (-Cj-A)/(B+Dj)/2.\n";
	
	// first row
	Cj = gsl_vector_get(C, 0);
	Dj = gsl_vector_get(D, 0);
	gsl_spmatrix_set(Coeff_Mat, 0,   0, (Dj+2.*A)/Dj);
	gsl_spmatrix_set(Coeff_Mat, 0, 0+1,     -2.*A/Dj);
	
	//out_h << " (Dj+A+B)/Dj     -2.*A/Dj\n";
	//out_h << setprecision(4) << (Dj+A)/Dj << " " 
	//	  << setprecision(4) <<  -2.*A/Dj << "\n";
	
	// interior
	for(int i = 1; i < M-1; ++i){
		Cj = gsl_vector_get(C, i);
		Dj = gsl_vector_get(D, i);
		
		gsl_spmatrix_set(Coeff_Mat, i, i-1,  (Cj-A)/Dj/2.);
		gsl_spmatrix_set(Coeff_Mat, i,   i,     (Dj+A)/Dj);
		gsl_spmatrix_set(Coeff_Mat, i, i+1, (-Cj-A)/Dj/2.);
		
		//out_h 	<< setprecision(4) <<  (Cj-A)/Dj/2. << " " 
				//<< setprecision(4) <<   (Dj+A)/Dj << " " 
				//<< setprecision(4) << (-Cj-A)/Dj/2. << "\n";
	}
	
	// last row
	Cj = gsl_vector_get(C, M-1);
	Dj = gsl_vector_get(D, M-1);
	gsl_spmatrix_set(Coeff_Mat, M-1, M-2, 0.);
	gsl_spmatrix_set(Coeff_Mat, M-1, M-1, 0.);
	

	// "convert to compressed column format"
	Coeff_Mat_Compressed = gsl_spmatrix_ccs(Coeff_Mat);
	
	
	//*******************************************************************************
	// THE SOLVER: GMRES
	
	// initialize x and y
		// define current rhs vector
		gsl_vector_memcpy(y, Tin);
		
		// initial guess for solution vector
		gsl_vector_memcpy(x, Tin);
	
	// Output the Inital Temperature Vector To Solution File
	for(int i = 0; i < M; ++i){
		out_d << setprecision(8) << gsl_vector_get(Tin, i) << " ";
	}
	out_d << "\n";
	
	// initiate and allocate memory for solver
	const double tol = 1.0e-8;	// "solution relative tolerance"
	const size_t max_iter = 10; // "maximum iterations"
	const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
	gsl_splinalg_itersolve *work = gsl_splinalg_itersolve_alloc(T, M, 100);
	size_t iter = 0;
	int status; // will return the status of the iterative solver employed by GSL
	
	//*******************************************************************************
	// MARCH ALONG SOLVING
	for(int i = 1; i < N; ++i){
		
		iter = 0;
		
		do {// iterate and return status
			status = gsl_splinalg_itersolve_iterate
				(Coeff_Mat_Compressed, y, tol, x, work);
		}
		while(status ==  GSL_CONTINUE && ++iter < max_iter);
		
		// Check to see if the GSL routine is iterating a lot
		if(iter == max_iter){ cout << "i = " << i << " Max iterations reached \n"; }
		
		// output the solution vector 
		for(int j = 0; j < M; ++j){
			out_d << setprecision(8) << gsl_vector_get(x, j) << " ";
		}
		out_d << "\n";
		
		// set y to x, as x is the new rhs vector for the next step forward
		// don't change x, the current step is the initial guess for the next step
		gsl_vector_memcpy(y, x);
	}
	
	//*******************************************************************************
	// clean up
	
	// free workspace allocated by gsl
	gsl_splinalg_itersolve_free(work);
	gsl_spmatrix_free(Coeff_Mat);
	gsl_spmatrix_free(Coeff_Mat_Compressed);
	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_vector_free(C);
	gsl_vector_free(D);
	gsl_vector_free(r);
	gsl_vector_free(u);
	
	// Close Output File Stream
	out_d.close();
	out_h.close();
	
	return 0;
}//end of main