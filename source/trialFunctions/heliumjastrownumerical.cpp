#include "heliumjastrownumerical.h"

HeliumJastrowNumerical::HeliumJastrowNumerical(VMCSolver *solver)
{
    simpleFlag = false;
    m_outfileName = "HeliumJastrowNumerical";

    solver->setCharge(2);
    solver->setNParticles(2);
    solver->setAlpha(1.843);
    solver->setBeta(0.34);

}

HeliumJastrowNumerical::~HeliumJastrowNumerical()
{

}

double HeliumJastrowNumerical::waveFunction(const mat &r, VMCSolver *solver)
{
    double nParticles = solver->getNParticles();
    double nDimensions = solver->getNDimensions();
    double alpha = solver->getAlpha();
    double beta = solver->getBeta();



    //double r12;
    vec rpos(nParticles);
    for(int i = 0; i < nParticles; i++) {
        double rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        rpos[i] = sqrt(rSingleParticle);
    }
    // assuming 2 particles
    //   (ta fra elektron-elektron pot.)
    //r12 = abs(rpos[0] - rpos[1]);
    double r12 = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = i + 1; j < nParticles; j++) {
            r12 = 0;
            for(int k = 0; k < nDimensions; k++) {
                r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
            }
        }
    }

    return exp(-accu(rpos) * alpha) * exp(r12 / (2.0*(1 + beta * r12))) ;
}

double HeliumJastrowNumerical::localEnergy(const mat &r, VMCSolver *solver)
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

