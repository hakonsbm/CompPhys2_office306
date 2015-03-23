#ifndef LIB_H
#define LIB_H
#include <armadillo>
using namespace arma;
#define   ZERO       1.0E-10
// Random number generator
double ran2(long *);

// function for gaussian random numbers
double GaussianDeviate(long *);

// decomposition of matrix
void ludcmp(mat &, int, int *, double *);

#endif // LIB_H


