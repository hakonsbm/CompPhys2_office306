#ifndef NEON_H
#define NEON_H

#include "trialfunction.h"
#include "vmcsolver.h"

#include <armadillo>

using namespace arma;

class VMCSolver;

class Neon: public TrialFunction
{

public:
    Neon(VMCSolver *solver);
    virtual double waveFunction(const mat &r, VMCSolver*  solver);
    virtual double localEnergy(const mat &r, VMCSolver *solver );
    virtual double lnDerivativeWaveFunction(const mat &r, VMCSolver *solver);
    virtual double lnSecondDerivativeWaveFunction(const mat &r, VMCSolver *solver);

private:
    double psi1s(double ri, double alpha); // Ansatz functions
    double psi2s(double ri, double alpha); // only for testing.
    double phi(const mat &r, double alpha, int i, int j, VMCSolver *solver);
    double SlaterDeterminant(const mat &r,double alpha, VMCSolver *solver);
    mat detUp;
    mat detDown;
};

#endif // NEON_H
