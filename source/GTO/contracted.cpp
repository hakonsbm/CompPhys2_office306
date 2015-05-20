// calculates primitive using weight and exponent
// calculates contracted GTO as sum of primitives, each with their own normalization

#include <armadillo>
#include <math.h>
#include "contracted.h"
#include "basisbank.h"

using namespace arma;
using namespace std;

Contracted::Contracted():
	// class variables here
	contracted(0)
{
}

// adds together the primitives to make contracted GTO
void Contracted::add_primitive(double alpha, double w, int i, int j, int k, const vec pos)
{
	double x = pos[0];
	double y = pos[1];
	double z = pos[2];
	contracted += normalization(alpha, i, j, k) * w * pow(x,i)*pow(y,j)*pow(z,k) * exp(-alpha*(x*x+y*y+z*z));
}

// normalization of GTO product
double Contracted::normalization(double alpha, int i, int j, int k)
{
	return pow(2*alpha/M_PI, 3/((double) 4)) * sqrt( pow(8*alpha,i+j+k) * 
			(fac(i)*fac(j)*fac(k)) / (fac(2*i)*fac(2*j)*fac(2*k)) )
}

// factorial function (because libraries are lacking)
int Contracted::fac(int n)
{
	int fac = 1;
	for (int i = 2; i <= n; ++i) fac = fac * i;
	return fac;
}
