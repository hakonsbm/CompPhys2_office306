#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <armadillo>
#include "trialFunctions/trialfunction.h"
#include "trialFunctions/heliumsimplenumerical.h"
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
    void calculateOptimalSteplength();
    void setTrialFunction(TrialFunction *trialFunction) { m_trialFunction = trialFunction; }
    void setAlpha(double alpha) {m_alpha = alpha;}
    void setBeta(double beta) {m_beta = beta;}
    void setSteplength(double SL) {stepLength = SL;}

    TrialFunction *trialFunction(){return m_trialFunction;}

    int getNDimensions() {return nDimensions; }
    double getAlpha() {return m_alpha; }
    double getBeta() {return m_beta; }
    int getCharge() {return charge; }
    double getH()   {return h;}
    double getH2()  {return h2;}
    int getNParticles() {return nParticles; }






private:
    TrialFunction *m_trialFunction;

    double waveFunction(const mat &r);
    double localEnergy(const mat &r);
    double QuantumForce(const mat &r, mat QForce);

    int nParticles;


    int nDimensions;
    int charge;
    double stepLength;
    double step_min;
    double m_alpha;
    double m_beta;
    double h;
    double h2;
    long idum;
    int nCycles;
    double D; // diffusion constant
    // double timestep; // timestep for gaussian deviate (using steplength)
    double GreensFunction;


    mat rOld;
    mat rNew;
    mat QForceOld;
    mat QForceNew;


};

#endif // VMCSOLVER_H
