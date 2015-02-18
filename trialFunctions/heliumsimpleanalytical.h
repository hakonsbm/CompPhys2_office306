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
    HeliumSimpleAnalytical();
    virtual double waveFunction(const mat &r, VMCSolver*  solver);
    virtual double localEnergy(const mat &r, VMCSolver *solver );

};

#endif // HELIUMSIMPLEANALYTICAL_H
