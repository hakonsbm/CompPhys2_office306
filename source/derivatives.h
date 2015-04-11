#ifndef DERIVATIVES_H
#define DERIVATIVES_H


#include <armadillo>

using namespace arma;

class VMCSolver;

class Derivatives
{
public:
    Derivatives();
    ~Derivatives();
    double numericalDoubleDerivative(const mat &r, VMCSolver *solver);

    double analyticalSimpleDoubleDerivative(const mat &r, VMCSolver *solver);

};

#endif // DERIVATIVES_H
