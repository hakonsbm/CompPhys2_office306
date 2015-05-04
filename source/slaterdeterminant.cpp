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
    //i is th particleTag and j is the type of wavefunction
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
    //i is th particleTag and j is the type of wavefunction
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
    detUpOld = zeros<mat>(nHalf, nHalf);
    detDownOld = zeros<mat>(nHalf, nHalf);


    for (int k = 0; k <  nHalf; ++k)
    {
        for (i = 0; i < nHalf; ++i)
        {
            // for detUp
            detUpOld(i,k) =  phi(r, alpha, i, k, solver);

            // for detDownOld
            detDownOld(i,k) =  phi(r, alpha, i + nHalf, k, solver);
        }
    }

    //If we are solving it analytically we also need the inverse of the slater matrix
    if(solver->trialFunction()->m_analytical)
    {
        detUpInverseOld = zeros<mat>(nHalf, nHalf);
        detDownInverseOld = zeros<mat>(nHalf, nHalf);

        detUpInverseOld = detUpOld.i();
        detDownInverseOld = detDownOld.i();
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
    mat tempdetDownOld = zeros<mat>(Nhalf, Nhalf);

    // fill matrix detUp and detDownOld
//    updateSlaterMatrices(r, solver);


    tempDetUp = detUpOld;
    tempdetDownOld = detDownOld;

    // decompose A (phi matrix) to B & C
    /*
     * End up with
     *     (c00 c01 c02 c03)
     * A = (b10 c11 c12 c13)
     *     (b20 b21 c22 c23)
     *     (b30 b31 b32 c33)
     */

    ludcmp(tempDetUp, Nhalf, indx, &d1);
    ludcmp(tempdetDownOld, Nhalf, indx, &d2);

    // compute SD as c00*c11*..*cnn
    SD = 1;
    for (i = 0; i < Nhalf; ++i)
    {
        SD *= tempDetUp(i, i)*tempdetDownOld(i, i);
    }
    // return SD
    return d1*d2*SD;
}


mat SlaterDeterminant::gradientSlaterDeterminant(const mat &r , VMCSolver *solver)
{

    //Not functinoining properly, needs to be several vectors since we need a seperate tracking of "force" on each particle,
    // 3 coordinates per particle

    int nParticles= solver->getNParticles();
    int nDimensions = solver->getNDimensions();
    mat gradient = zeros (nParticles, nDimensions);

    int nHalf = nParticles/2;

//    updateSlaterMatrices(r,solver);

    //Calculating the sum of the particles derivatives
    for(int i = 0; i < nHalf; i ++) //Sums over the particles
    {
        for(int j = 0; j < nHalf; j++)
        {
            gradient.row(i) +=  gradientPhi(r,i,j,solver).t()*detUpInverseOld(j,i);
            gradient.row(i + nHalf) += gradientPhi(r, i + nHalf, j, solver).t() * detDownInverseOld(j,i);
        }
    }

    return gradient;

}
 
double SlaterDeterminant::laplacianSlaterDeterminant(const mat &r, VMCSolver *solver)
{
    double derivative = 0;
    int nParticles= solver->getNParticles();
    int nHalf = nParticles/2;

    updateSlaterMatrices(r,solver);

    //Calculating the sum of the particles derivatives
    for(int i = 0; i < nHalf; i ++) //Sums over the particles
    {
        for(int j = 0; j < nHalf; j++)
        {
            derivative += laplacianPhi(r, i, j, solver) * detUpInverseOld(j,i);

            derivative += laplacianPhi(r, i + nHalf, j, solver) * detDownInverseOld(j,i);
        }
    }

    return derivative;

}
