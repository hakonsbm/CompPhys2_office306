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

private:
    double psi1s(double ri, double alpha); // Ansatz functions
    double psi2s(double ri, double alpha); // only for testing.
    double phi(double ri, double alpha, int M, const vec &ri_vec);
    double SlaterDeterminant(const mat &r, int nParticles, int nDimensions, double alpha);
    mat detUp;
    mat detDown;
};

#endif // NEON_H
