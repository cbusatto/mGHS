/* File: utils_rgamma3p.h */

#ifndef utils_rgamma3p_H
#define utils_rgamma3p_H

// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// Parabolic Cylinder D

double rg3p (const double& a, const double& b, const int& c);
double rg3p_approx (const double& a, const double& b, const int& c);
double rg3p_c1 (const double& a, const double& b);
double rg3p_c1_fast (const double& a, const double& b);

#endif

