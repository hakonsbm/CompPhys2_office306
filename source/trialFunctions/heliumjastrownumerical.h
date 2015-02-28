#ifndef HELIUMJASTROWNUMERICAL_H
#define HELIUMJASTROWNUMERICAL_H

#include "trialfunction.h"
#include "../vmcsolver.h"

#include <armadillo>

using namespace arma;

class VMCSolver;

class HeliumJastrowNumerical : public TrialFunction
{
public:
    HeliumJastrowNumerical(VMCSolver *solver);
    ~HeliumJastrowNumerical();

    double waveFunction(const mat &r, VMCSolver *solver);
    double localEnergy(const mat &r, VMCSolver *solver);

};

#endif // HELIUMJASTROWNUMERICAL_H
