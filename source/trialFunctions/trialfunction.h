#ifndef TRIALFUNCTION_H
#define TRIALFUNCTION_H

#include <armadillo>

using namespace arma;

class VMCSolver;

class TrialFunction
{
public:
    TrialFunction();
    virtual double waveFunction(const mat &r, VMCSolver *solver ) = 0 ;
    virtual double localEnergy(const mat &r, VMCSolver *solver ) = 0;

    bool simpleFlag;

};

#endif // TRIALFUNCTION_H
