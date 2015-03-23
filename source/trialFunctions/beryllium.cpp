#include "beryllium.h"
#include "trialfunction.h"
#include "vmcsolver.h"
#include "lib.h"

#include <iostream>

using namespace std;

Beryllium::Beryllium(VMCSolver *solver)
{
    simpleFlag = false;
    m_outfileName = "Beryllium";

    solver->setCharge(4);
    solver->setNParticles(4);
    solver->setAlpha(4.0);
    solver->setBeta(0.31);
}

double Beryllium::waveFunction(const mat &r, VMCSolver *solver)
{
    double rSingleParticle, alpha, beta, wf, SD, product, rij, a;
    int spin_count = 0;
    vec spins(6);
    spins(0) = 1./2.; spins(5) = 1./2.;
    for(int s = 1; s<=4; s++) spins(s) = 1./4.;
    product = 1.0;
    alpha = solver -> getAlpha();
    beta = solver -> getBeta();
    //vec argument(solver->getNParticles());
    for(int i = 0; i < solver->getNParticles(); i++) {
        //argument[i] = 0.0;
        rSingleParticle = 0;
        /*
        for(int j = 0; j < solver->getNDimensions(); j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        argument[i] = sqrt(rSingleParticle);
        */
        for(int j = i + 1; j < solver->getNParticles(); j++) {
            rij = 0;
            for(int k = 0; k < solver->getNDimensions(); k++) {
                rij += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
            }
            rij = sqrt(rij);
            a = spins(spin_count);
            product = product * exp(a*rij/(1+beta*rij));
            spin_count++;
        }
    }
    /*
    wf = (psi1s(argument[0], alpha)*psi2s(argument[1], alpha)
        -psi1s(argument[1], alpha)*psi2s(argument[0], alpha))*
        (psi1s(argument[2], alpha)*psi2s(argument[3], alpha)
        -psi1s(argument[3], alpha)*psi2s(argument[2], alpha));
    */
    SD = SlaterDeterminant(r, solver->getNParticles(), solver->getNDimensions(), alpha);
    //cout << "wf / SD: " << wf << " / " << SD << endl; // check if we get the expected value
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

double Beryllium::phi(double ri, double alpha, int M)
{
    // returns an ansatz based on matrix row, M, in the SD
    if (M == 0)
    {
        return exp(-alpha*ri); // 1s
    }
    else if (M == 1)
    {
        return (1-alpha*ri/2.0)*exp(-alpha*ri/2.0); // 2s
    }
    else if (M>=2 && M<=4)
    {
        return alpha*ri*exp(-alpha*ri/2.0); // 2p
    }
}




double Beryllium::SlaterDeterminant(const mat &r, int nParticles,
                                    int nDimensions, double alpha)
{
    int i, j, Nhalf, M, *indx;
    double d1, d2, SD, ri;
    Nhalf = nParticles/2;
    indx = new int [Nhalf];
    detUp = zeros<mat>(Nhalf, Nhalf);
    detDown = zeros<mat>(Nhalf, Nhalf);
    // fill matrix detUp and detDown
    for (M = 0; M <  Nhalf; ++M)
    {
        for (i = 0; i < Nhalf; ++i)
        {
            // for detUp
            ri = 0;
            for (j = 0; j < nDimensions; ++j) ri += r(i,j)*r(i,j);
            ri = sqrt(ri);
            detUp(i,M) =  phi(ri, alpha, M);

            // for detDown
            ri = 0;
            for (j = 0; j < nDimensions; ++j) ri += r(i + Nhalf,j)*r(i + Nhalf,j);
            ri = sqrt(ri);
            detDown(i,M) =  phi(ri, alpha, M);
        }
    }
    // decompose A (phi matrix) to B & C
    /*
     * End up with
     *     (c00 c01 c02 c03)
     * A = (b10 c11 c12 c13)
     *     (b20 b21 c22 c23)
     *     (b30 b31 b32 c33)
     */
    ludcmp(detUp, Nhalf, indx, &d1);
    ludcmp(detDown, Nhalf, indx, &d2);
    // compute SD as c00*c11*..*cnn
    SD = 1;
    for (i = 0; i < Nhalf; ++i)
    {
        SD *= detUp(i, i)*detDown(i, i);
    }
    // return SD
    return d1*d2*SD;
}
