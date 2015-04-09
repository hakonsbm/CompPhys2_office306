#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <armadillo>
#include "trialFunctions/trialfunction.h"
#include "trialFunctions/heliumsimplenumerical.h"
#include "trialFunctions/heliumsimpleanalytical.h"
#include "trialFunctions/heliumjastrownumerical.h"
#include "trialFunctions/heliumjastrowanalytical.h"
#include "derivatives.h"

using namespace arma;

class TrialFunction; class Derivatives;

class VMCSolver
{
public:
    VMCSolver();

    void runMasterIntegration();
    void runMonteCarloIntegration();
    void runMonteCarloIntegrationIS();

    void calculateOptimalSteplength();
    void QuantumForce(const mat &r, mat &QForce);
    void setTrialFunction(TrialFunction *trialFunction) { m_trialFunction = trialFunction; }
    void setAlpha(double alpha) {m_alpha = alpha;}
    void setBeta(double beta) {m_beta = beta;}
    void setCharge(int C) {charge = C;}
    void setNParticles(int P) {nParticles = P;}
    void setStepLength(double inStepLength) {stepLength = inStepLength;}
    void setCycles(double cycles) {nCycles = cycles; }

    void initiateDerivatives(Derivatives *derivatives) {m_derivatives = derivatives; }
    Derivatives *derivatives(){return m_derivatives;}

    TrialFunction *trialFunction(){return m_trialFunction;}

    int getNParticles() {return nParticles; }
    int getNDimensions() {return nDimensions; }
    double getAlpha() {return m_alpha; }
    double getBeta() {return m_beta; }
    int getCharge() {return charge; }
    double getH()   {return h;}
    double getH2()  {return h2;}
    double getStepLength() {return stepLength;}
    int getMy_Rank() {return my_rank;}
    void switchbBlockSampling(bool onOff) { m_blockSampling = onOff;}


    double getEnergyVar() {return m_energyVar;}
    double getEnergy() {return m_energy;}
    void storeEnergy(double energy) {m_energy = energy;}
    void storeVariance(double variance) {m_energyVar = variance;}

//    void mpiArguments( int nargs, char* args[]){ m_nargs = nargs; m_args = args; }

private:
    TrialFunction *m_trialFunction;
    Derivatives *m_derivatives;

    double waveFunction(const mat &r);
    double localEnergy(const mat &r);

    bool importanceSampling;    //When this flag is true it uses importance sampling instead of regular sampling
    bool m_blockSampling;


    int nDimensions;
    int charge;
    double stepLength;
    double step_min;
    double m_alpha;
    double m_beta;
    int nParticles;
    double h;
    double h2;
    long idum;
    int nCycles;

    //Variables for the MPI implementation
    int my_rank;
//    int m_nargs;
//    char** m_args;

    //Each thread should store data as common accessible doubles, so they can be merged into the the master thread
    double m_energyVar;
    double m_energy;
    double m_energySquared;
    double m_averageR12;
    double m_ratio;
    double m_moves;

    double D; // diffusion constant
    // double timestep; // timestep for gaussian deviate (using steplength)
    double GreensFunction;


    mat rOld;
    mat rNew;
    mat QForceOld;
    mat QForceNew;


};

#endif // VMCSOLVER_H
