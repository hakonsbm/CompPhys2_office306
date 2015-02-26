#include "beryllium.h"
#include "trialfunction.h"
#include "vmcsolver.h"


#include <iostream>

using namespace std;

Beryllium::Beryllium()
{
    simpleFlag = false;
    m_outfileName = "Beryllium";
}

double Beryllium::waveFunction(const mat &r, VMCSolver *solver)
{

    double rSingleParticle, alpha, beta, wf, product, rij;
    product = 1.0;
    alpha = solver -> getAlpha();
    beta = solver -> getBeta();
    vec argument(solver->getNParticles());
    for(int i = 0; i < solver->getNParticles(); i++) {
        argument[i] = 0.0;
        rSingleParticle = 0;
        
        for(int j = 0; j < solver->getNDimensions(); j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        argument[i] = sqrt(rSingleParticle);
        for(int j = i + 1; j < solver->getNParticles(); j++) {
            rij = 0;
            for(int k = 0; k < solver->getNDimensions(); k++) {
                rij += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
            }
            product = product * exp(rij/(2.0*(1+beta*rij)));
        }
    }

    wf = (psi1s(argument[0], alpha)*psi2s(argument[1], alpha)
        -psi1s(argument[1], alpha)*psi2s(argument[0], alpha))*
        (psi1s(argument[2], alpha)*psi2s(argument[3], alpha)
        -psi1s(argument[3], alpha)*psi2s(argument[2], alpha));

    return wf*product;
}

double Beryllium::localEnergy(const mat &r, VMCSolver *solver)
{
    //Grabbing all the necessary constants stored in the solver
    double nParticles = solver->getNParticles();
    double nDimensions = solver->getNDimensions();
    double charge = solver->getCharge();
    double h = solver->getH();
    double h2 = solver->getH2();


    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);

    rPlus = rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    double waveFunctionCurrent = waveFunction(r, solver);



    // Kinetic energy

    double kineticEnergy = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = waveFunction(rMinus, solver);
            waveFunctionPlus = waveFunction(rPlus, solver);
            kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }
    kineticEnergy = 0.5 * h2 * kineticEnergy / waveFunctionCurrent;

    // Potential energy
    double potentialEnergy = 0;
    double rSingleParticle = 0;
    for(int i = 0; i < nParticles; i++) {
        rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j)*r(i,j);
        }
        potentialEnergy -= charge / sqrt(rSingleParticle);
    }
    // Contribution from electron-electron potential
    double r12 = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = i + 1; j < nParticles; j++) {
            r12 = 0;
            for(int k = 0; k < nDimensions; k++) {
                r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
            }
            potentialEnergy += 1 / sqrt(r12);
        }
    }

    return kineticEnergy + potentialEnergy;
}

double Beryllium::psi1s(double ri, double alpha)
{
    return exp(-alpha*ri);
}


double Beryllium::psi2s(double ri, double alpha)
{
    return (1-alpha*ri/2.0)*exp(-alpha*ri/2.0);
}
