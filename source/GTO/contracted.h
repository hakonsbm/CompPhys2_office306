#ifndef Contracted_H
#define Contracted_H

#include <armadillo>

using namespace arma;

class Contracted{
	
public:
    Contracted();
    void add_primitive(double alpha, double w, int i, int j, int k, const vec pos);
    double get_orb_1s() {return orb_1s;}
    double get_orb_2s() {return orb_2s;}
    double get_orb_pX() {return orb_pX;}
    double get_orb_pY() {return orb_pY;}
    double get_orb_pZ() {return orb_pZ;}

    void contract_orb_1s();
    void contract_orb_2s();
    void contract_orb_pX();
    void contract_orb_pY();
    void contract_orb_pZ();

    void contract_orb();
    void contract_orbK();

private:
    double contracted;
    double orb_1s, orb_2s, orb_pX, orb_pY, orb_pZ;
    double orb, K;
	double normalization(double alpha, int i, int j, int k);
    int fac(int n);
};
#endif // Contracted_H
