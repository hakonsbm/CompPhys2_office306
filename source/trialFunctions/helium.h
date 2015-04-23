#ifndef HELIUM_H
#define HELIUM_H


#include "trialfunction.h"
#include "../vmcsolver.h"

#include <armadillo>

using namespace arma;

class VMCSolver;

class Helium : public TrialFunction
{

public:
    Helium(VMCSolver* solver);
    virtual double waveFunction(const mat &r, VMCSolver*  solver);
    virtual double localEnergy(const mat &r, VMCSolver *solver );


};

#endif // HELIUM_H
