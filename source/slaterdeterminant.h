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

    double calculateDeterminant(const mat &r,double alpha, VMCSolver *solver);
    double phi(const mat &r, double alpha, int i, int j, VMCSolver *solver);

    mat detUp;
    mat detDown;
};

#endif // SLATERDETERMINANT_H
