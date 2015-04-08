#ifndef DERIVATIVES_H
#define DERIVATIVES_H

class VMCSolver;

class Derivatives
{
public:
    Derivatives();
    ~Derivatives();
    double numericalDerivative(VMCSolver *solver);
};

#endif // DERIVATIVES_H
