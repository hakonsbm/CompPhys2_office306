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

};

#endif // BERYLLIUM_H
