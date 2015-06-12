#ifndef HELIUMJASTROWANALYTICAL_H
#define HELIUMJASTROWANALYTICAL_H

#include "trialfunction.h"
#include "../vmcsolver.h"

#include <armadillo>

using namespace arma;

class VMCSolver;

class HeliumJastrowAnalytical : public TrialFunction
{
public:
    HeliumJastrowAnalytical(VMCSolver *solver);
    ~HeliumJastrowAnalytical();

    virtual double waveFunction(const mat &r, VMCSolver *solver);
    virtual double localEnergy(const mat &r, VMCSolver *solver);
    virtual double lnDerivativeWaveFunction(const mat &r, VMCSolver *solver );
    virtual double lnSecondDerivativeWaveFunction(const mat &r, VMCSolver *solver );

};

#endif // HELIUMJASTROWANALYTICAL_H
