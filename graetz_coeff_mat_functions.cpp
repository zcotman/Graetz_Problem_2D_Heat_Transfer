// graetz_coeff_mat_functions.cpp
//
// For use with: graetz_problem_conduct.cpp
//
// Zachary Cotman
//
// 11/21/2018
//
// These functions are used to build the coefficient matricies
//	used in the GSL sparse-matrix solver GMRES

#include <cmath>
using namespace std;

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>


//***********************************************************************************
// Function Declarations
void set_coeff_interior
(const int M, const int A, const int B,
const gsl_vector *C, const gsl_vector *D, gsl_spmatrix *Coeff_Interior);

void set_coeff_final
(const int M, const int A, const int B,
const gsl_vector *C, const gsl_vector *D, gsl_spmatrix *Coeff_Final);

//***********************************************************************************
// Function Definitions
void set_coeff_interior
(const int M, const int A, const int B,
const gsl_vector *C, const gsl_vector *D, gsl_spmatrix *Coeff_Interior)
{
	double Cj = 0, Dj = 0; Cj = Dj;
	
		// top row
		Cj = gsl_vector_get(C, 0);
		Dj = gsl_vector_get(D, 0);
		gsl_spmatrix_set(Coeff_Interior, 0, 0, 4.*A + 2.*B + Dj);
		gsl_spmatrix_set(Coeff_Interior, 0, 1,       -4.*A);
		
		// interior rows
		for(int j = 1; j < M-1; ++j){
			Cj = gsl_vector_get(C, j);
			Dj = gsl_vector_get(D, j);
		
			gsl_spmatrix_set(Coeff_Interior, j, j-1,     Cj-A);
			gsl_spmatrix_set(Coeff_Interior, j,   j, 2.*(A+B)+Dj);
			gsl_spmatrix_set(Coeff_Interior, j, j+1,    -Cj-A);
		}
		
		// bottom row
		Cj = gsl_vector_get(C, M-1);
		Dj = gsl_vector_get(D, M-1);
		gsl_spmatrix_set(Coeff_Interior, M-1, M-2, 0.);
		gsl_spmatrix_set(Coeff_Interior, M-1, M-1, 0.);
		
		// last row is all 0s
}
//***********************************************************************************
void set_coeff_final
(const int M, const int A, const int B,
const gsl_vector *C, const gsl_vector *D, gsl_spmatrix *Coeff_Final)
{
	double Cj = 0, Dj = 0;
	
	// first row
	Cj = gsl_vector_get(C, 0);
	Dj = gsl_vector_get(D, 0);
	gsl_spmatrix_set(Coeff_Final, 0, 0, 4.*A + 2.*B + Dj);
	gsl_spmatrix_set(Coeff_Final, 0, 1,      -4.*A);	
	
	// interior rows
	for(int j = 1; j < M-1; ++j){
		Cj = gsl_vector_get(C, j);
		Dj = gsl_vector_get(D, j);
		
		gsl_spmatrix_set(Coeff_Final, j, j+1,        Cj-A);
		gsl_spmatrix_set(Coeff_Final, j,   j, 2.*(A+B)+Dj);
		gsl_spmatrix_set(Coeff_Final, j, j-1,       -Cj-A);
	}
	
	// bottom row
		Cj = gsl_vector_get(C, M-1);
		Dj = gsl_vector_get(D, M-1);
		gsl_spmatrix_set(Coeff_Final, M-1, M-2, 0.);
		gsl_spmatrix_set(Coeff_Final, M-1, M-1, 0.);
	
	// last row is all 0s
}