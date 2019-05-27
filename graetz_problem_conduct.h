// graetz_problem_conduct.h
//
// Zachary Cotman
//
// 11/21/2018
//
// For linking the external functions to graetz_problem_conduct.cpp

#ifndef GRAETZ_PROBLEM_CONDUCT_H
#define GRAETZ_PROBLEM_CONDUCT_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>

// ***************************************************************************
// From graetz_residual_functions.cpp
extern double 
res_interior
	(const int N, const int M, const double A, const double B,
	const gsl_vector * C, const gsl_vector * D, const gsl_matrix * T);

extern double 
res_bottom
	(const int N, const int M, const double A, const double B,
	const gsl_vector * C, const gsl_vector * D, const gsl_matrix * T);

extern double 
res_lastline
	(const int N, const int M, const double A, const double B,
	const gsl_vector * C, const gsl_vector * D, const gsl_matrix * T);

extern double 
res_corner
	(const int N, const int M, const double A, const double B,
	const gsl_vector * C, const gsl_vector * D, const gsl_matrix * T);
	
// ***************************************************************************
// From graetz_coeff_mat_functions.cpp
extern void
set_coeff_interior
	(const int M, const int A, const int B,
	const gsl_vector *C, const gsl_vector *D, gsl_spmatrix *Coeff_Interior);

extern void
set_coeff_final
	(const int M, const int A, const int B,
	const gsl_vector *C, const gsl_vector *D, gsl_spmatrix *Coeff_Final);

// ***************************************************************************
// From graetz_xy_functions.cpp
extern void 
set_initial_guess
	(const int M, const int i, gsl_vector *x, const gsl_matrix *T);

extern void 
update_T
	(const int M, const int i, const gsl_vector *x, gsl_matrix *T);

extern void 
set_rhs_interior
	(const int M, const int i, const double B,
	const gsl_vector *D, const gsl_matrix *T, gsl_vector *y);

extern void 
set_rhs_final
	(const int M, const int i, const double B,
	const gsl_vector *D, const gsl_matrix *T, gsl_vector *y);

#endif