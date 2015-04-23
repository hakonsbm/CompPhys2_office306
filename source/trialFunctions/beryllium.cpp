#include "beryllium.h"
#include "trialfunction.h"
#include "vmcsolver.h"
#include "lib.h"

#include <iostream>

using namespace std;

Beryllium::Beryllium(VMCSolver *solver)
{
    simpleFlag = false;
    m_analytical = false;
    m_outfileName = "Beryllium";

    solver->setCharge(4);
    solver->setNParticles(4);
    solver->setAlpha(4.0);
    solver->setBeta(0.31);

    //Giving the particles in Beryllium it's spin, the first half up and the second part down
    // up = 0 and down = 1
    spin << 0 << 0 << 1 << 1;



}

double Beryllium::waveFunction(const mat &r, VMCSolver *solver)
{
    double rSingleParticle, alpha, beta, wf, SD, product, rij, a;
    product = 1.0;
    alpha = solver -> getAlpha();
    beta = solver -> getBeta();
    //vec argument(solver->getNParticles());
    if(solver->getElectronInteration())
    {
        for(int i = 0; i < solver->getNParticles(); i++) {
            rSingleParticle = 0;
            for(int j = i + 1; j < solver->getNParticles(); j++) {
                rij = 0;
                for(int k = 0; k < solver->getNDimensions(); k++) {
                    rij += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
                }
                rij = sqrt(rij);
                a = spinFactor(i,j);
                product = product * exp(a*rij/(1+beta*rij));
            }
        }
    }

    SD = solver->determinant()->calculateDeterminant(r,alpha,solver); //SlaterDeterminant(r, alpha, solver);

    return SD*product;
}

double Beryllium::localEnergy(const mat &r, VMCSolver *solver)
{
    //Grabbing all the necessary constants stored in the solver
    double nParticles = solver->getNParticles();
    double nDimensions = solver->getNDimensions();
    double charge = solver->getCharge();
    vec gradientSlater, gradientJastrow;

    // Kinetic energy

    double kineticEnergy = 0;

    if(m_analytical)
    {   //Calculates the kinetic energy as the ratios of -1/2* ( d²/dx²|D| /|D| + 2 * (d/dx |D|/|D|)*d/dx Psi_C/Psi_C + d²/dx² Psi_C /Psi_C )
        kineticEnergy += solver->determinant()->laplacianSlaterDeterminant(r, solver);
//        gradientJastrow = solver->derivatives()->analyticalCorrelationDerivative(r,solver);
        gradientSlater = solver->determinant()->gradientSlaterDeterminant(r,solver);
//        kineticEnergy += 2*(dot(solver->determinant()->gradientSlaterDeterminant(r,solver), solver->derivatives()->analyticalCorrelationDerivative(r,solver));
        kineticEnergy += solver->derivatives()->analyticalCorrelationDoubleDerivative(r,solver);
        kineticEnergy *= -1./2.;
    }
    else
        kineticEnergy = solver->derivatives()->numericalDoubleDerivative(r, solver) / (2.*waveFunction(r, solver));

//    cout << endl<< "KineticEnergy by the numerical: " << solver->derivatives()->numericalDoubleDerivative(r, solver) / (2.*waveFunction(r, solver)) << endl;

//    cout << " KineticEnergy by the analytical: " << solver->determinant()->laplacianSlaterDeterminant(r, solver)/(-2.) << endl << endl;


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
    if(solver->getElectronInteration())
    {
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
    }

    return kineticEnergy + potentialEnergy;
}

double Beryllium::spinFactor(int i, int j)     //The corrolation factor a in the Jastrow factor 1/2 if opposite spin or 1/4 if same
{
    if(spin(i) == spin(j))
        return 1./4.;
    else
        return 1./2.;
}


//double Beryllium::psi1s(double ri, double alpha)
//{
//    return exp(-alpha*ri);
//}


//double Beryllium::psi2s(double ri, double alpha)
//{
//    return (1-alpha*ri/2.0)*exp(-alpha*ri/2.0);
//}

//double Beryllium::phi(const mat &r, double alpha, int i, int j, VMCSolver *solver)
//{
//    // returns an ansatz based on matrix row, M, in the SD
//    int nDimensions = solver->getNDimensions();
//    double ri = 0;
//    for (int k = 0; k < nDimensions; ++k) ri += r(i,k)*r(i,k);
//    ri = sqrt(ri);

//    if (j == 0)
//    {
//        return exp(-alpha*ri); // 1s
//    }
//    else if (j == 1)
//    {
//        return (1-alpha*ri/2.0)*exp(-alpha*ri/2.0); // 2s
//    }
//    else if (j>=2 && j<=4)
//    {
//        int dimension = j-2;
//        return alpha*r(i,dimension)*exp(-alpha*ri/2.0); // 2p
//    }
//}



//double Beryllium::SlaterDeterminant(const mat &r,double alpha, VMCSolver *solver)
//{
//    int i, j, Nhalf, *indx;
//    double d1, d2, SD;
//    int nParticles= solver->getNParticles();
//    Nhalf = nParticles/2;
//    indx = new int [Nhalf];
//    detUp = zeros<mat>(Nhalf, Nhalf);
//    detDown = zeros<mat>(Nhalf, Nhalf);
//    // fill matrix detUp and detDown

//    for (int k = 0; k <  Nhalf; ++k)
//    {
//        for (i = 0; i < Nhalf; ++i)
//        {
//            // for detUp
//            detUp(i,k) =  phi(r, alpha, i, k, solver);

//            // for detDown
//            detDown(i,k) =  phi(r, alpha, i+Nhalf, k, solver);
//        }
//    }
//    // decompose A (phi matrix) to B & C
//    /*
//     * End up with
//     *     (c00 c01 c02 c03)
//     * A = (b10 c11 c12 c13)
//     *     (b20 b21 c22 c23)
//     *     (b30 b31 b32 c33)
//     */

//    ludcmp(detUp, Nhalf, indx, &d1);
//    ludcmp(detDown, Nhalf, indx, &d2);

//    // compute SD as c00*c11*..*cnn
//    SD = 1;
//    for (i = 0; i < Nhalf; ++i)
//    {
//        SD *= detUp(i, i)*detDown(i, i);
//    }
//    // return SD
//    return d1*d2*SD;
//}
