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
    void numericalGradient(mat &gradient, const mat &r, VMCSolver *solver);

    void analyticalGradient(mat &gradient, const mat &r, VMCSolver *solver);
    void analyticalLaplacianRatio(double &laplacianRatio, const mat &r, VMCSolver *solver);

//    double analyticalCorrelationDoubleDerivative(const mat &r, VMCSolver *solver);

    vec analyticalPsi1SDerivative(int particleTag, const mat &r, VMCSolver *solver);
    double analyticalPsi1SDoubleDerivative(int particleTag, const mat &r, VMCSolver *solver);
    vec analyticalPsi2SDerivative(int particleTag, const mat &r, VMCSolver *solver);
    double analyticalPsi2SDoubleDerivative(int particleTag, const mat &r, VMCSolver *solver);
    vec analyticalPsi2PDerivative(int particleTag, int dimension, const mat &r, VMCSolver *solver);
    double analyticalPsi2PDoubleDerivative(int particleTag, int dimension, const mat &r, VMCSolver *solver);

    mat  analyticalCorrelationGradient(const mat &r, VMCSolver *solver);
    double analyticalCorrelationDoubleDerivative( const mat &r, VMCSolver *solver);
    double analyticalCorrelationLaplacian( const mat &r, VMCSolver *solver);


    double fDerivative(int i, int j, const mat &r, VMCSolver *solver);
    double fDoubleDerivative(int i, int j, const mat &r, VMCSolver *solver);
};

#endif // DERIVATIVES_H
