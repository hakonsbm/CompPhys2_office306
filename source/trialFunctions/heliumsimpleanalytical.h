#ifndef HELIUMSIMPLEANALYTICAL_H
#define HELIUMSIMPLEANALYTICAL_H


#include "trialfunction.h"
#include "../vmcsolver.h"

#include <armadillo>

using namespace arma;

class VMCSolver;

class HeliumSimpleAnalytical : public TrialFunction
{

public:
    HeliumSimpleAnalytical(VMCSolver* solver);
    virtual double waveFunction(const mat &r, VMCSolver*  solver);
    virtual double localEnergy(const mat &r, VMCSolver *solver );
    virtual double lnDerivativeWaveFunction(const mat &r, VMCSolver *solver);
    virtual double lnSecondDerivativeWaveFunction(const mat &r, VMCSolver *solver);


};

#endif // HELIUMSIMPLEANALYTICAL_H
