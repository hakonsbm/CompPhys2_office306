#ifndef GTO_H
#define GTO_H

#include <armadillo>
#include "basisbank.h"

using namespace arma;

class GTO{
public:
    GTO();
    double GTO_phi(string sys, const mat &r, int i, int j);
    //void initBasis(basisbank *basis) {m_basis = basis;}
    //basisbank *basis() {return m_basis;}

private:
    //basisbank *m_basis;
	/*
	void Contracted_orbital_1s();
	void Contracted_orbital_2s();
	void Contracted_orbital_3s();
	void Contracted_orbital_1p();
	void Contracted_orbital_2p();
	*/
};
#endif // GTO_H
