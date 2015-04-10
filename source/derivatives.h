#ifndef DERIVATIVES_H
#define DERIVATIVES_H

class VMCSolver;

class Derivatives
{
public:
    Derivatives();
    ~Derivatives();
    double numericalDoubleDerivative(VMCSolver *solver);
};

#endif // DERIVATIVES_H
