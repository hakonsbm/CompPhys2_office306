#ifndef BERYLLIUMTWO_H
#define BERYLLIUMTWO_H

#include "trialfunction.h"
#include "vmcsolver.h"

#include <armadillo>

using namespace arma;

class VMCSolver;

class BerylliumTwo: public TrialFunction
{

public:
    BerylliumTwo(VMCSolver *solver);
    virtual double waveFunction(const mat &r, VMCSolver*  solver);
    virtual double localEnergy(const mat &r, VMCSolver *solver );
    virtual double lnDerivativeWaveFunction(const mat &r, VMCSolver *solver);
    virtual double lnSecondDerivativeWaveFunction(const mat &r, VMCSolver *solver);
//    double spinFactor(int i, int j);     //The corrolation factor a in the Jastrow factor 1/2 if opposite spin or 1/4 if same

private:

};



#endif // BERYLLIUMTWO_H
