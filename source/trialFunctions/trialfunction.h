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
    virtual double waveFunction(const mat &r, VMCSolver *solver ) = 0;
    virtual double localEnergy(const mat &r, VMCSolver *solver ) = 0;
    virtual double lnDerivativeWaveFunction(const mat &r, VMCSolver *solver ) = 0;
    virtual double lnSecondDerivativeWaveFunction(const mat &r, VMCSolver *solver ) = 0;
    void setAnalytical(bool onOff) {m_analytical = onOff; }
    double spinFactor(int i, int j);
    double getNucleusDistance() {return m_nucleusDistance;}
    void setNucleusDistance(double R) ;
    void calculateAlpha(VMCSolver *solver);  //This is only for the cases where alpha is directly calculatable, for example H_2 and Be_2, it calculates and sets alpha
    void setConjugate(bool onOff) { m_conjugateMethod = onOff; }
    bool getConjugate() {return m_conjugateMethod; }


    string m_outfileName;

    bool simpleFlag;
    bool m_analytical;
    bool m_conjugateMethod;
    bool m_zeroDistance;    //For use with the molecules, so 1/|R| is not included if R = 0

    ivec spin;

private:
    double m_nucleusDistance;


};

#endif // TRIALFUNCTION_H
