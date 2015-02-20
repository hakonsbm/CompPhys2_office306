#ifndef HELIUMSIMPLENUMERICALLY_H
#define HELIUMSIMPLENUMERICALLY_H

#include "trialfunction.h"
#include "../vmcsolver.h"

#include <armadillo>

using namespace arma;

class VMCSolver;

class HeliumSimpleNumerical : public TrialFunction
{

public:
    HeliumSimpleNumerical();
    virtual double waveFunction(const mat &r, VMCSolver*  solver);
    virtual double localEnergy(const mat &r, VMCSolver *solver );

};

#endif // HELIUMSIMPLENUMERICALLY_H
