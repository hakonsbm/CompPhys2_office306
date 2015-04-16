#ifndef SLATERDETERMINANT_H
#define SLATERDETERMINANT_H
#include <armadillo>
#include "vmcsolver.h"
#include "derivatives.h"


using namespace arma;

class VMCSolver;

class SlaterDeterminant
{
public:
    SlaterDeterminant();
    ~SlaterDeterminant();

    void updateSlaterMatrices(const mat &r, VMCSolver *solver);
    double calculateDeterminant(const mat &r,double alpha, VMCSolver *solver);
    double phi(const mat &r, double alpha, int i, int j, VMCSolver *solver);
    double SlaterDeterminant::determinantRatioUp(const mat &r, VMCSolver *solver, Derivatives *der);
    double SlaterDeterminant::determinantRatioDown(const mat &r, VMCSolver *solver, Derivatives *der);
    double SlaterDeterminant::determinantLaplacianRatioUp(const mat &r, VMCSolver *solver, Derivatives *der);
    double SlaterDeterminant::determinantLaplacianRatioDown(const mat &r, VMCSolver *solver, Derivatives *der);
    double laplacianPhi( const mat &r, double alpha, int i, int j, VMCSolver *solver);
    double gradientSlaterDeterminant(const mat &r , VMCSolver *solver);
    double laplacianSlaterDeterminant(const mat &r , VMCSolver *solver);

    mat detUp;
    mat detDown;
    mat detUpInverse;
    mat detDownInverse;
};

#endif // SLATERDETERMINANT_H
