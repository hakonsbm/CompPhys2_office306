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
    HeliumJastrowAnalytical();
    ~HeliumJastrowAnalytical();

    double waveFunction(const mat &r, VMCSolver *solver);
    double localEnergy(const mat &r, VMCSolver *solver);

};

#endif // HELIUMJASTROWANALYTICAL_H
