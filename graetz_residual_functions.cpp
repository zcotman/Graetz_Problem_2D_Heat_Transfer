// graetz_residual_functions.cpp
//
// For use with: graetz_problem_conduct.cpp
//
// Zachary Cotman
//
// 11/21/2018
//
// These functions are used to calculate the residuals of the
//	iterative solver


#include <fstream>
#include <cmath>
using namespace std; // meh

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>
//#include <gsl/gsl_splinalg.h>

//#include "graetz_problem_conduct.h"

//***********************************************************************************
// Function Prototypes 
double res_interior
(const int N, const int M, const double A, const double B,
const gsl_vector * C, const gsl_vector * D, const gsl_matrix * T);

double res_bottom
(const int N, const int M, const double A, const double B,
const gsl_vector * C, const gsl_vector * D, const gsl_matrix * T);

double res_lastline
(const int N, const int M, const double A, const double B,
const gsl_vector * C, const gsl_vector * D, const gsl_matrix * T);

double res_corner
(const int N, const int M, const double A, const double B,
const gsl_vector * C, const gsl_vector * D, const gsl_matrix * T);
					
//***********************************************************************************
// Function Definitions
double res_interior	(const int N, const int M, const double A, const double B,
					const gsl_vector * C, const gsl_vector * D, const gsl_matrix * T)
{
	double res = 0;
	double res_tot = 0;
	double Cj = 0;
	double Dj = 0;
	
	for(int i = 1; i < N-1; ++i){
		for(int j = 1; j < M-1; ++j){
				
			Cj = gsl_vector_get(C, j);
			Dj = gsl_vector_get(D, j);
	
			res = 	(Cj - A)       *gsl_matrix_get(T, j-1,   i) + 
					(2.*(A + B)+Dj)*gsl_matrix_get(T,   j,   i) + 
					(-Cj - A)      *gsl_matrix_get(T, j+1,   i) -
					(B + Dj)       *gsl_matrix_get(T,   j, i-1) - 
					B              *gsl_matrix_get(T,   j, i+1);
				
		res_tot += res*res;
		
		}
	}
	
	return res_tot;
}		

//***********************************************************************************
double res_bottom	(const int N, const int M, const double A, const double B,
					const gsl_vector * C, const gsl_vector * D, const gsl_matrix * T)
{
	double res = 0;
	double res_tot = 0;
	int j = M; j = 0;
	double Cj = gsl_vector_get(C, j);
	double Dj = gsl_vector_get(D, j);
	j = 0;
	
	for(int i = 1; i < N-1; ++i){
		res = 	(4.*A + 2.*B + Dj)*gsl_matrix_get(T,   j,   i) +
				-4.*A             *gsl_matrix_get(T, j+1,   i) -
				(B + Dj)          *gsl_matrix_get(T,   j, i-1) - 
				B                 *gsl_matrix_get(T,   j, i+1);
		res_tot += res*res;
	}
	
	return res_tot;
}

//***********************************************************************************
double res_lastline	(const int N, const int M, const double A, const double B,
					const gsl_vector * C, const gsl_vector * D, const gsl_matrix * T)
{
	double res = 0;
	double res_tot = 0;
	double Cj = 0;
	double Dj = 0;
	int i = N-1;
	
	for(int j = 1; j < M-1; ++j){
		
		Cj = gsl_vector_get(C, j);
		Dj = gsl_vector_get(D, j);
	
		res = 	(Cj - A)       *gsl_matrix_get(T, j-1,   i) + 
				(2.*(A + B)+Dj)*gsl_matrix_get(T,   j,   i) + 
				(-Cj - A)      *gsl_matrix_get(T, j+1,   i) -
				2.*(B + Dj)    *gsl_matrix_get(T,   j, i-1);
				
		res_tot += res*res;
	}
	
	return res_tot;
}

//***********************************************************************************
double res_corner	(const int N, const int M, const double A, const double B,
					const gsl_vector * C, const gsl_vector * D, const gsl_matrix * T)
{
	double res = 0;
	double res_tot = 0;
	int i = N-1;
	int j = M-M; j = 0;
	double Cj = gsl_vector_get(C, j);
	double Dj = gsl_vector_get(D, j);
	
	res = 	(4.*A + 2.*B + Dj)*gsl_matrix_get(T,   j,   i) + 
			-4.*A             *gsl_matrix_get(T, j+1,   i) -
			2.*(B + Dj)       *gsl_matrix_get(T,   j, i-1);
			
	res_tot = res*res;
	
	return res_tot;
}