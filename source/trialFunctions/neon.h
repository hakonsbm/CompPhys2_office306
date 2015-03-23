#ifndef NEON_H
#define NEON_H

#include "trialfunction.h"
#include "../vmcsolver.h"
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
    double psi1s(double ri, double alpha);
    double psi2s(double ri, double alpha);
    double psi2p(double ri, double alpha);
};

#endif // NEON_H
