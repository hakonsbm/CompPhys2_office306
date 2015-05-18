#ifndef HYDROGENTWO_H
#define HYDROGENTWO_H

#include "trialfunction.h"
#include "../vmcsolver.h"

#include <armadillo>

using namespace arma;

class VMCSolver;

class HydrogenTwo : public TrialFunction
{
public:

    ~HydrogenTwo();

    HydrogenTwo(VMCSolver* solver);
    virtual double waveFunction(const mat &r, VMCSolver*  solver);
    virtual double localEnergy(const mat &r, VMCSolver *solver );


};

#endif // HYDROGENTWO_H
