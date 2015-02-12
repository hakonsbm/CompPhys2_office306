#ifndef TRIALFUNCTION_H
#define TRIALFUNCTION_H

#include <armadillo>

using namespace arma;

class VMCSolver;

class TrialFunction
{
public:
    TrialFunction();
    double m_localEnergy;
//    void testFunction( VMCSolver *solver );
    virtual double waveFunction(const mat &r, VMCSolver *solver ) = 0 ;

};

#endif // TRIALFUNCTION_H
