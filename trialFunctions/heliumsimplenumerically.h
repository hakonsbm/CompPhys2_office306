#ifndef HELIUMSIMPLENUMERICALLY_H
#define HELIUMSIMPLENUMERICALLY_H

#include "trialfunction.h"
#include "../vmcsolver.h"

#include <armadillo>

using namespace arma;

class VMMSolver;

class HeliumSimpleNumerically : public TrialFunction
{

public:
    HeliumSimpleNumerically();
    double localEnergy;
    virtual double waveFunction(const mat &r, VMCSolver*  solver);

};

#endif // HELIUMSIMPLENUMERICALLY_H
