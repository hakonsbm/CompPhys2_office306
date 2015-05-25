// calculates primitive using weight and exponent
// calculates contracted GTO as sum of primitives, each with their own normalization

#include <armadillo>
#include <math.h>
#include <ostream>
#include "contracted.h"
#include "basisbank.h"

using namespace arma;
using namespace std;

Contracted::Contracted():
	// class variables here
	contracted(0),
    orb(0),
	orb_1s(0),
	orb_2s(0),
    orb_pX(0),
    orb_pY(0),
    orb_pZ(0),
    K(0.5)
{
}

// adds together the primitives to make contracted GTO
void Contracted::add_primitive(double alpha, double w, int i, int j, int k, const vec pos)
{
	double x = pos(0);
	double y = pos(1);
	double z = pos(2);
    //cout << "pos: " << pos << endl;
    //cout << "alpha=" << alpha << "  w="<< w << endl;
	contracted += normalization(alpha, i, j, k) * w * pow(x,i)*pow(y,j)*pow(z,k) * exp(-alpha*(x*x+y*y+z*z));
    //cout << "primitive= " << normalization(alpha, i, j, k) * w * pow(x,i)*pow(y,j)*pow(z,k) * exp(-alpha*(x*x+y*y+z*z)) << "   ";
    //cout << "norm= " << normalization(alpha, i, j, k)<< "   ";
    //cout << "cont = " << contracted << endl;
}

// normalization of GTO product
double Contracted::normalization(double alpha, int i, int j, int k)
{
    /*
    cout << "alpha = " << alpha << "  ";
    cout << "2alph/pi = " << pow(2*alpha/M_PI, 3/((double) 4)) << "  ";
    cout << "8alph = " << pow(8*alpha,i+j+k) << "  ";
    cout << "fac = " << fac(i)*fac(j)*fac(k) << "  ";
    cout << "fac2 = " << (fac(2*i)*fac(2*j)*fac(2*k)) << "  ";
    cout << "N = " << pow(2*alpha/M_PI, 3/((double) 4)) * sqrt( pow(8*alpha,i+j+k) * (fac(i)*fac(j)*fac(k)) / (fac(2*i)*fac(2*j)*fac(2*k)) ) << endl;
    */
    return pow(2*alpha/M_PI, 3/((double) 4)) * sqrt( pow(8*alpha,i+j+k) *
            (fac(i)*fac(j)*fac(k)) / (fac(2*i)*fac(2*j)*fac(2*k)) );
}

// factorial function (because libraries are lacking)
int Contracted::fac(int n)
{
	int fac = 1;
	for (int i = 2; i <= n; ++i) fac = fac * i;
	return fac;
}

// contraction of orbitals
void Contracted::contract_orb_1s()
{
    orb_1s += orb + contracted*K;
    //cout << "orb_1s = " << contracted << endl;
	contracted = 0;
}

void Contracted::contract_orb_2s()
{
    orb_2s = contracted;
    //cout << "orb_2s= " << contracted << endl;
	contracted = 0;
}

void Contracted::contract_orb_pX()
{
    orb_pX += contracted*K;
	contracted = 0;
}

void Contracted::contract_orb_pY()
{
    orb_pY += contracted*K;
	contracted = 0;
}

void Contracted::contract_orb_pZ()
{
    orb_pZ += contracted*K;
    contracted = 0;
}

void Contracted::contract_orb()
{
    orb += contracted;
    contracted = 0;
}

void Contracted::contract_orbK()
{
    orb += contracted*K;
    contracted = 0;
}


