#ifndef TRIALFUNCTION_H
#define TRIALFUNCTION_H

#include <armadillo>

using namespace arma;
using namespace std;

class VMCSolver;

class TrialFunction
{
public:
    TrialFunction();
    virtual double waveFunction(const mat &r, VMCSolver *solver ) = 0 ;
    virtual double localEnergy(const mat &r, VMCSolver *solver ) = 0;
    void setAnalytical(bool onOff) {m_analytical = onOff; }

    string m_outfileName;

    bool simpleFlag;
    bool m_analytical;



};

#endif // TRIALFUNCTION_H
