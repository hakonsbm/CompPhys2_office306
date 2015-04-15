#ifndef SLATERDETERMINANT_H
#define SLATERDETERMINANT_H
#include <armadillo>


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
    double laplacianPhi( const mat &r, double alpha, int i, int j, VMCSolver *solver);
    double gradientSlaterDeterminant(const mat &r , VMCSolver *solver);
    double laplacianSlaterDeterminant(const mat &r , VMCSolver *solver);

    mat detUp;
    mat detDown;
    mat detUpInverse;
    mat detDownInverse;
};

#endif // SLATERDETERMINANT_H
