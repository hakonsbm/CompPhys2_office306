#include "slaterdeterminant.h"
#include "vmcsolver.h"
#include "lib.h"

SlaterDeterminant::SlaterDeterminant()
{

}

SlaterDeterminant::~SlaterDeterminant()
{

}

double SlaterDeterminant::phi(const mat &r, double alpha, int i, int j, VMCSolver *solver)
{
    // returns an ansatz based on matrix row, M, in the SD
    int nDimensions = solver->getNDimensions();
    double ri = 0;
    for (int k = 0; k < nDimensions; ++k) ri += r(i,k)*r(i,k);
    ri = sqrt(ri);

    if (j == 0)
    {
        return exp(-alpha*ri); // 1s
    }
    else if (j == 1)
    {
        return (1-alpha*ri/2.0)*exp(-alpha*ri/2.0); // 2s
    }
    else if (j>=2 && j<=4)
    {
        int dimension = j-2;
        return alpha*r(i,dimension)*exp(-alpha*ri/2.0); // 2p
    }
}


vec SlaterDeterminant::gradientPhi(const mat &r, int i, int j, VMCSolver *solver)
{
    vec derivative;

    if (j == 0)
    {
        derivative = solver->derivatives()->analyticalPsi1SDerivative(i,r,solver);
        return derivative;  // d²/dx² 1s
    }
    else if (j == 1)
    {
        derivative = solver->derivatives()->analyticalPsi2SDerivative(i,r,solver);

        return derivative; // d²/dx²  2s
    }
    else if (j>=2 && j<=4)
    {
        int dimension = j-2;
        derivative = solver->derivatives()->analyticalPsi2PDerivative(i,dimension,r,solver);

        return derivative; // d²/dx²  2p
    }
}

double SlaterDeterminant::laplacianPhi(const mat &r, int i, int j, VMCSolver *solver)
{
    double derivative;

    if (j == 0)
    {
        derivative = solver->derivatives()->analyticalPsi1SDoubleDerivative(i,r,solver);
        return derivative;  // d²/dx² 1s
    }
    else if (j == 1)
    {
        derivative = solver->derivatives()->analyticalPsi2SDoubleDerivative(i,r,solver);

        return derivative; // d²/dx²  2s
    }
    else if (j>=2 && j<=4)
    {
        int dimension = j-2;
        derivative = solver->derivatives()->analyticalPsi2PDoubleDerivative(i,dimension,r,solver);

        return derivative; // d²/dx²  2p
    }
}

void SlaterDeterminant::updateSlaterMatrices(const mat &r, VMCSolver *solver)
{
    int i;
    int nHalf= solver->getNParticles()/2;
    double alpha = solver->getAlpha();
    detUp = zeros<mat>(nHalf, nHalf);
    detDown = zeros<mat>(nHalf, nHalf);



    for (int k = 0; k <  nHalf; ++k)
    {
        for (i = 0; i < nHalf; ++i)
        {
            // for detUp
            detUp(i,k) =  phi(r, alpha, i, k, solver);

            // for detDown
            detDown(i,k) =  phi(r, alpha, i + nHalf, k, solver);
        }
    }

    //If we are solving it analytically we also need the inverse of the slater matrix
    if(solver->trialFunction()->m_analytical)
    {
        detUpInverse = zeros<mat>(nHalf, nHalf);
        detDownInverse = zeros<mat>(nHalf, nHalf);

        detUpInverse = detUp.i();
        detDownInverse = detDown.i();
    }

}

double SlaterDeterminant::calculateDeterminant(const mat &r,double alpha, VMCSolver *solver)
{
    int i, j, Nhalf, *indx;
    double d1, d2, SD;
    int nParticles= solver->getNParticles();
    Nhalf = nParticles/2;
    indx = new int [Nhalf];
    mat tempDetUp = zeros<mat>(Nhalf, Nhalf);
    mat tempDetDown = zeros<mat>(Nhalf, Nhalf);

    // fill matrix detUp and detDown
    updateSlaterMatrices(r, solver);

    tempDetUp = detUp;
    tempDetDown = detDown;

    // decompose A (phi matrix) to B & C
    /*
     * End up with
     *     (c00 c01 c02 c03)
     * A = (b10 c11 c12 c13)
     *     (b20 b21 c22 c23)
     *     (b30 b31 b32 c33)
     */

    ludcmp(tempDetUp, Nhalf, indx, &d1);
    ludcmp(tempDetDown, Nhalf, indx, &d2);

    // compute SD as c00*c11*..*cnn
    SD = 1;
    for (i = 0; i < Nhalf; ++i)
    {
        SD *= tempDetUp(i, i)*tempDetDown(i, i);
    }
    // return SD
    return d1*d2*SD;
}

//double SlaterDeterminant::determinantRatioUp(const mat &r, VMCSolver *solver, Derivatives *der)
//{
//    double determinantRatio = 0;
//    int nHalf = solver->getNParticles() / 2;
//    mat inverse = detUp.i();
//    for(int i = 0; i < nHalf; i++) {
//        for(int j = 0; j < nHalf; j++) {
//            if(j == 0) {
//                determinantRatio += der->analyticalPsi1SDerivative(i, &r, *solver) * inverse(j,i);
//            }
//            else if(j == 1) {
//                determinantRatio += der->analyticalPsi2SDerivative(i, &r, *solver) * inverse(j,i);
//            }
//            else if((j >= 2) && (nHalf > 2)) {
//                determinantRatio += der->analyticalPsi2PDerivative(i, &r, *solver) * inverse(j,i);
//            }
//        }
//    }
//    return determinantRatio / solver->getMHR();
//}

//double SlaterDeterminant::determinantRatioDown(const mat &r, VMCSolver *solver, Derivatives *der)
//{
//    double determinantRatio = 0;
//    int nHalf = solver->getNParticles() / 2;
//    mat inverse = detDown.i();
//    for(int i = 0; i < nHalf; i++) {
//        for(int j = 0; j < nHalf; j++) {
//            if(j == 0) {
//                determinantRatio += der->analyticalPsi1SDerivative(i + nHalf, &r, *solver) * inverse(j,i);
//            }
//            else if(j == 1) {
//                determinantRatio += der->analyticalPsi2SDerivative(i + nHalf, &r, *solver) * inverse(j,i);
//            }
//            else if((j >= 2) && (nHalf > 2)) {
//                determinantRatio += der->analyticalPsi2PDerivative(i + nHalfi, &r, *solver) * inverse(j,i);
//            }
//        }
//    }
//    return determinantRatio / solver->getMHR();
//}

vec SlaterDeterminant::gradientSlaterDeterminant(const mat &r , VMCSolver *solver)
{
    vec derivative = zeros (3);
    int nParticles= solver->getNParticles();
    int nHalf = nParticles/2;

//    updateSlaterMatrices(r,solver);

<<<<<<< HEAD
    //Calculating the sum of the particles derivatives
    for(int i = 0; i < nHalf; i ++) //Sums over the particles
    {
        for(int j = 0; j < nHalf; j++)
        {
            derivative += gradientPhi(r, i, j, solver) * detUpInverse(j,i);

            derivative += gradientPhi(r, i + nHalf, j, solver) * detDownInverse(j,i);
        }
    }

    return derivative;
=======
vec SlaterDeterminant::gradientSlaterDeterminant(const mat &r , VMCSolver *solver)
{
       //Shall calculate (d/dx |D|)/|D|

    return zeros(3);
>>>>>>> be6af1ebcd312b68aeb0c74e650ec16d881f998f
}

double SlaterDeterminant::laplacianSlaterDeterminant(const mat &r, VMCSolver *solver)
{
    double derivative = 0;
    int nParticles= solver->getNParticles();
    int nHalf = nParticles/2;

//    updateSlaterMatrices(r,solver);

    //Calculating the sum of the particles derivatives
    for(int i = 0; i < nHalf; i ++) //Sums over the particles
    {
        for(int j = 0; j < nHalf; j++)
        {
            derivative += laplacianPhi(r, i, j, solver) * detUpInverse(j,i);

            derivative += laplacianPhi(r, i + nHalf, j, solver) * detDownInverse(j,i);
        }
    }

    return derivative;

}
