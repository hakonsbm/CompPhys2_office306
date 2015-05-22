// Calculates contracted GTO orbitals as sum of primitives, each with their own normalization.
// Primitives are calculated using weight and exponent.
// Works only in 3 dimensions.

#include <armadillo>
#include <math.h>
#include "gto.h"
#include "basisbank.h"

using namespace arma;
using namespace std;

// calls basisbank which sends weights and exponents to be contracted
GTO::GTO(string sys, const mat &r, int i)
	// class variables here
{
    initBasis(new basisbank);
    vec pos(3);
	pos(0)=r(i,0); pos(1)=r(i,1); pos(2)=r(i,2);

    if (sys == "he") basis()->add_3_21G_he(pos); // for helium atom
    else if (sys == "ne") basis()->add_3_21G_ne(pos); // for neon atom
    else if (sys == "h") basis()->add_3_21G_h(pos); // for hydrogen atom
    else if (sys == "be") basis()->add_3_21G_be(pos); // for beryllium atom
}

// returns orbital based on matrix row, M in the SD
double GTO::GTO_phi(int j)
{
    if (j == 0) return basis()->get_orb_1s(); // 1s
    else if (j == 1) return basis()->get_orb_2s(); // 2s
    else if (j == 2) return basis()->get_orb_pX(); // 2p
    else if (j == 3) return basis()->get_orb_pY();
    else if (j == 4) return basis()->get_orb_pZ();

}

/*
void GTO::Contracted_orbital_1s()
{
	return basis->get_orb_1s();
}

void GTO::Contracted_orbital_2s()
{
	return basis->get_orb_2s();
}

void GTO::Contracted_orbital_3s()
{
	return basis->get_orb_3s();
}

void GTO::Contracted_orbital_1p()
{
	return basis->get_orb_1p();
}

void GTO::Contracted_orbital_2p()
{
	return basis->get_orb_2p();
}

*/
