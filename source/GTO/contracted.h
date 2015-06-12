#ifndef Contracted_H
#define Contracted_H

#include <armadillo>

using namespace arma;

class Contracted{
	
public:
    Contracted();
    void add_primitive(double alpha, double w, int i, int j, int k, const vec pos);
    double get_orb() {return orb;}
    void contract_orb(double coeff);

private:
    double contracted;
    double orb;
	double normalization(double alpha, int i, int j, int k);
    int fac(int n);
};
#endif // Contracted_H
