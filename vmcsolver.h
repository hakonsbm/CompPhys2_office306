#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <armadillo>
#include "trialFunctions/trialfunction.h"
#include "trialFunctions/heliumsimplenumerically.h"

using namespace arma;

class TrialFunction; class HeliumSimpleNumerically;

class VMCSolver
{
public:
    VMCSolver();

    void runMonteCarloIntegration();

    void setTrialFunction(TrialFunction *trialFunction) { m_trialFunction = trialFunction; }
    TrialFunction *trialFunction(){return m_trialFunction;}

    int getNParticleshere() {return nParticles; }
    int getNDimensions() {return nDimensions; }
    double getAlpha() {return alpha; }
    double getBeta() {return beta; }
    int getCharge() {return charge; }


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
