// calculates primitive using weight and exponent
// calculates contracted GTO as sum of primitives, each with their own normalization

#include <armadillo>
#include <math.h>
#include "gto.h"
#include "basisbank.h"

using namespace arma;
using namespace std;

// calls basisbank which sends weights and exponents to be contracted
GTO::GTO(char sys, const vec pos):
	// class variables here
{
	basisbank basis = new basisbank();
	if (sys == "he") basis->add_3_21G_he(const vec pos); // for helium atom
	if (sys == "ne") basis->add_3_21G_ne(const vec pos); // for neon atom
	if (sys == "h") basis->add_3_21G_h(const vec pos); // for hydrogen atom
}

GTO::Contracted_orbital()
{
	return basis->get_contracted();
}
