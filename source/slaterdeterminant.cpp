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

double SlaterDeterminant::calculateDeterminant(const mat &r,double alpha, VMCSolver *solver)
{
    int i, j, Nhalf, *indx;
    double d1, d2, SD;
    int nParticles= solver->getNParticles();
    Nhalf = nParticles/2;
    indx = new int [Nhalf];
    detUp = zeros<mat>(Nhalf, Nhalf);
    detDown = zeros<mat>(Nhalf, Nhalf);
    // fill matrix detUp and detDown

    for (int k = 0; k <  Nhalf; ++k)
    {
        for (i = 0; i < Nhalf; ++i)
        {
            // for detUp

            detUp(i,k) =  phi(r, alpha, i, k, solver);

            // for detDown
            detDown(i,k) =  phi(r, alpha, i+Nhalf, k, solver);
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
