// graetz_xy_functions.cpp
//
// For use with: graetz_problem_conduct.cpp
//
// Zachary Cotman
//
// 11/21/2018
//
// These functions are used to build the right hand side(rhs), y, and
//	initial guess, x, vectors used in the GSL sparse-matrix solver GMRES

#include <cmath>
using namespace std; // meh

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>

//***********************************************************************************
// Function Declarations
void set_initial_guess
(const int i, gsl_vector *x, const gsl_matrix *T);

void update_T
(const int i, const gsl_vector *x, gsl_matrix *T);

void set_rhs_interior
(const int M, const int i, const double B,
const gsl_vector *D, gsl_vector *y, const gsl_matrix *T);

void set_rhs_final
(const int M, const int i, const double B,
const gsl_vector *D, gsl_vector *y, const gsl_matrix *T);

//***********************************************************************************
// Function Definitions
void set_initial_guess
(const int M, const int i, gsl_vector *x, const gsl_matrix *T)
{
	for(int j = 0; j < M; ++j){
		gsl_vector_set(x, j, gsl_matrix_get(T, j, i));
	}	
}

//***********************************************************************************
void update_T
(const int M, const int i, const gsl_vector *x, gsl_matrix *T)
{
	for(int j = 0; j < M; ++j){
		gsl_matrix_set(T, j, i, gsl_vector_get(x, j));
	}
}

//***********************************************************************************
void set_rhs_interior
(const int M, const int i, const double B, 
const gsl_vector *D, const gsl_matrix *T, gsl_vector *y)
{
	double Dj = 0;
	
	for(int j = 0; j < M; ++j){
		Dj = gsl_vector_get(D, j);
		gsl_vector_set(y, j,  (B+Dj)*gsl_matrix_get(T, j, i-1)
							     + B*gsl_matrix_get(T, j, i+1));
	}
}

//***********************************************************************************
void set_rhs_final
(const int M, const int i, const double B,
const gsl_vector *D, const gsl_matrix *T, gsl_vector *y)
{
	double Dj = 0;
	
	for(int j = 0; j < M; ++j){
		Dj = gsl_vector_get(D, j);
		gsl_vector_set(y, j, 2.*(B+Dj)*gsl_matrix_get(T, j, i-1));
	}
}