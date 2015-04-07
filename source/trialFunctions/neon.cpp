#include "neon.h"
#include "trialfunction.h"
#include "vmcsolver.h"
#include "libm.h"
#include <iostream>
#include <cmath>

using namespace std;

Neon::Neon(VMCSolver *solver)
{
    simpleFlag = false;
    m_outfileName = "Neon";

    solver->setCharge(10);
    solver->setNParticles(10);
    solver->setAlpha(4.0);//We don't know what should be here
    solver->setBeta(0.31);//or here either
}


double Neon::waveFunction(const mat &r, VMCSolver *solver)
{
    double rSingleParticle, alpha, beta, wf, product, rij, a;
    double slater[solver->getNParticles()][solver->getNParticles()];
    int trash[solver->getNParticles()];
    double moretrash;
    int spin_count = 0;
    vec spins(33);//What's this...?
    spins(0) = 1./2.; spins(33) = 1./2.;
    for(int s = 1; s<=32; s++) spins(s) = 1./4.;
    product = 1.0;
    wf = 1.0;
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
            a = spins(spin_count);
            product = product * exp(a*rij/(1+beta*rij));
            spin_count++;
        }
    }
    for(int i = 0; i < solver->getNParticles(); i++) {
        for(int j = 0; j < solver->getNDimensions(); j++) {
            if((i == 0) || (i == 1)) {
                slater[i][j] = psi1s(argument[j], alpha);}
            else if((i == 2) || (i == 3)) {
                slater[i][j] = psi2s(argument[j], alpha);}
            else if((i > 3) && (i < 10)) {
                slater[i][j] = psi2p(argument[j], alpha);}
        }
    }

//    ludcmp(slater, solver->getNParticles(), trash, moretrash);

    for(int i = 0; i < solver->getNParticles(); i++) {
        wf *= slater[i][i];
    }
    wf /= 720*sqrt(7);

    return wf*product;
}

double Neon::localEnergy(const mat &r, VMCSolver *solver)
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

double Neon::psi1s(double ri, double alpha)
{
    return exp(-alpha*ri);
}


double Neon::psi2s(double ri, double alpha)
{
    return (1-alpha*ri/2.0)*exp(-alpha*ri/2.0);
}

double Neon::psi2p(double ri, double alpha)
{
    return alpha*ri*exp(-alpha*ri/2.0);
}
