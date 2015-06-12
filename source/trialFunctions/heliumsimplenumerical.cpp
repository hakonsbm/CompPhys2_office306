#include "heliumsimplenumerical.h"
#include "trialfunction.h"
#include "vmcsolver.h"


#include <iostream>

using namespace std;

HeliumSimpleNumerical::HeliumSimpleNumerical(VMCSolver *solver)
{
    simpleFlag = true;
    m_outfileName = "HeliumSimpleNumerical";

    solver->setCharge(2);
    solver->setNParticles(2);
    solver->setAlpha(1.65);
    solver->setBeta(0);

}

double HeliumSimpleNumerical::waveFunction(const mat &r, VMCSolver *solver)
{

    //Calculates the wavefunction Psi = prod_i exp(-alpha *r_i)

    vec rpos(solver->getNParticles());
    for(int i = 0; i < solver->getNParticles(); i++) {
        double rSingleParticle = 0;
        for(int j = 0; j < solver->getNDimensions(); j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        rpos[i] = sqrt(rSingleParticle);
    }

    return exp(-accu(rpos) * solver->getAlpha());
}

double HeliumSimpleNumerical::localEnergy(const mat &r, VMCSolver *solver)
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

    kineticEnergy = solver->derivatives()->numericalDoubleDerivative(r, solver) / (2.*waveFunctionCurrent);



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

//     Contribution from electron-electron potential
    if(solver->getElectronInteration())
    {
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
    }

    return kineticEnergy + potentialEnergy;
}

double HeliumSimpleNumerical::lnDerivativeWaveFunction(const mat &r, VMCSolver *solver)
{
    return 0;
}

double HeliumSimpleNumerical::lnSecondDerivativeWaveFunction(const mat &r, VMCSolver *solver)
{
    return 0;
}
