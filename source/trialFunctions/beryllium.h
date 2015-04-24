#ifndef BERYLLIUM_H
#define BERYLLIUM_H

#include "trialfunction.h"
#include "vmcsolver.h"

#include <armadillo>

using namespace arma;

class VMCSolver;

class Beryllium: public TrialFunction
{

public:
    Beryllium(VMCSolver *solver);
    virtual double waveFunction(const mat &r, VMCSolver*  solver);
    virtual double localEnergy(const mat &r, VMCSolver *solver );
    double spinFactor(int i, int j);     //The corrolation factor a in the Jastrow factor 1/2 if opposite spin or 1/4 if same

private:
//    ivec spin;
//    double psi1s(double ri, double alpha); // Ansatz functions
//    double psi2s(double ri, double alpha); // only for testing.
//    double phi(const mat &r, double alpha, int i, int j, VMCSolver *solver);
//    double SlaterDeterminant(const mat &r,double alpha, VMCSolver *solver);
//    mat detUp;
//    mat detDown;
};

#endif // BERYLLIUM_H
