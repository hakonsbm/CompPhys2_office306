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

private:
    double psi1s(double ri, double alpha); // Ansatz functions
    double psi2s(double ri, double alpha); // only for testing.
    double phi(double ri, double alpha, int M);
    double SlaterDeterminant(const mat &r, int nParticles, int nDimensions, double alpha);
    mat detUp;
    mat detDown;
};

#endif // BERYLLIUM_H
