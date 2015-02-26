#ifndef BERYLLIUM_H
#define BERYLLIUM_H

#include "trialfunction.h"
#include "../vmcsolver.h"

#include <armadillo>

using namespace arma;

class VMCSolver;

class Beryllium: public TrialFunction
{

public:
    Beryllium();
    virtual double waveFunction(const mat &r, VMCSolver*  solver);
    virtual double localEnergy(const mat &r, VMCSolver *solver );

private:
    double psi1s(double ri, double alpha);
    double psi2s(double ri, double alpha);

};

#endif // BERYLLIUM_H
