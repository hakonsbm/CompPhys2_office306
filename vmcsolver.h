#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <armadillo>
#include "trialFunctions/trialfunction.h"
#include "trialFunctions/heliumsimplenumerically.h"
#include "trialFunctions/heliumsimpleanalytical.h"
#include "trialFunctions/heliumjastrownumerical.h"
#include "trialFunctions/heliumjastrowanalytical.h"

using namespace arma;

class TrialFunction; class HeliumSimpleNumerically;

class VMCSolver
{
public:
    VMCSolver();

    void runMonteCarloIntegration();

    void setTrialFunction(TrialFunction *trialFunction) { m_trialFunction = trialFunction; }
    TrialFunction *trialFunction(){return m_trialFunction;}

    int getNParticles() {return nParticles; }
    int getNDimensions() {return nDimensions; }
    double getAlpha() {return alpha; }
    double getBeta() {return beta; }
    int getCharge() {return charge; }
    double getH()   {return h;}
    double getH2()  {return h2;}


private:
    TrialFunction *m_trialFunction;

    double waveFunction(const mat &r);
    double localEnergy(const mat &r);

    int nDimensions;
    int charge;
    double stepLength;
    int nParticles;

    double h;
    double h2;

    long idum;

    double alpha;
    double beta;

    int nCycles;

    mat rOld;
    mat rNew;


};

#endif // VMCSOLVER_H
