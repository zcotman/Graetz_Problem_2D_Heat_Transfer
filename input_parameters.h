// input_parameters.h
//
// Zachary Cotman
//
// Nov 22, 2018
//
// Edit these parameters in order to change the simulation mesh
//
// SMALLER VALUES FOR dx AND dr MIGHT OVERLOAD RAM! see below
 
#ifndef INPUT_PARAMETERS_H
#define INPUT_PARAMETERS_H

#include <cmath>
const double PI = atan(1.)*4.;

const double alpha  = (PI*0.6/2./4182)*(PI*0.6/2./4182);
//const double alpha   = 0.1;

const int my_max_iter = 50000;

const double radius = 1.;
const double length = 0.31; // fully developed @ ~ 0.1
const double dx     = 0.0001; // >= 0.0001	
const double dr     = 0.02; // >= 0.0005	
const int N         = round(length/dx)+1;
const int M         = round(radius/dr)+1;
#endif