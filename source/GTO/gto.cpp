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
GTO::GTO()
	// class variables here
{
}

// returns orbital based on matrix row, M in the SD
double GTO::GTO_phi(string sys, const mat &r, int i, int j)
{
    //initBasis(basisbank);
    basisbank basis;
    vec pos(3);

    pos(0)=r(i,0); pos(1)=r(i,1); pos(2)=r(i,2);

    if (sys == "He"){
        vec coeffs(2);
        coeffs(0) = 0.4579;
        coeffs(1) = 0.6573;
        basis.add_3_21G_he(pos, coeffs); // for helium atom
        //cout <<  coeffs << endl;
        return basis.get_orb();
    }

    else if (sys == "Be"){
        vec coeffs(9);
        if (j == 0) {
            coeffs(0) = -9.9281e-01;
            coeffs(1) = -7.6425e-02;
            coeffs(2) = 2.8727e-02;
            coeffs(3) = 1.2898e-16;
            coeffs(4) = -2.3257e-19;
            coeffs(5) = 5.6097e-19;
            coeffs(6) = 1.2016e-16;
            coeffs(7) = -4.6874e-19;
            coeffs(8) = 1.1319e-18;

        }
        else if (j == 1){
            coeffs(0) = -2.1571e-01;
            coeffs(1) = 2.2934e-01;
            coeffs(2) = 8.2235e-01;
            coeffs(3) = 5.1721e-16;
            coeffs(4) = 4.5670e-18;
            coeffs(5) = -1.1040e-17;
            coeffs(6) = 8.5306e-16;
            coeffs(7) = 7.0721e-18;
            coeffs(8) = -1.7060e-17;
        }
        else{
            cout << "ERROR! WRONG J IN GTO";
            exit(1);
        }

        basis.add_3_21G_be(pos, coeffs); // for beryllium atom
        return basis.get_orb();
    }
    else if (sys == "Ne"){
        vec coeffs(9);

        if (j == 0) {
            coeffs(0) = -9.807700e-01;
            coeffs(1) = -9.371400e-02;
            coeffs(2) = 2.286300e-02;
            coeffs(3) = -9.951900e-19;
            coeffs(4) = -1.212500e-18;
            coeffs(5) = -4.180000e-19;
            coeffs(6) = -1.669600e-19;
            coeffs(7) = 1.212500e-18;
            coeffs(8) = 3.877900e-19;
        }
        else if (j == 1) {
            coeffs(0) = -2.606200e-01;
            coeffs(1) = 2.585800e-01;
            coeffs(2) = 8.161900e-01;
            coeffs(3) = -5.618600e-18;
            coeffs(4) = -2.861500e-16;
            coeffs(5) = 4.619900e-17;
            coeffs(6) = -4.240500e-18;
            coeffs(7) = -2.942600e-16;
            coeffs(8) = 5.051900e-17;
        }
        else if (j == 2) {
            coeffs(0) = 1.159600e-16;
            coeffs(1) = -2.010600e-16;
            coeffs(2) = -3.236100e-16;
            coeffs(3) = 2.715500e-02;
            coeffs(4) = -5.620700e-01;
            coeffs(5) = 9.113900e-03;
            coeffs(6) = 2.889000e-02;
            coeffs(7) = -5.979700e-01;
            coeffs(8) = 9.695900e-03;
        }
        else if (j == 3) {
            coeffs(0) = -8.371600e-18;
            coeffs(1) = -9.717300e-17;
            coeffs(2) = 1.323700e-16;
            coeffs(3) = -4.032000e-01;
            coeffs(4) = -2.583300e-02;
            coeffs(5) = -3.918000e-01;
            coeffs(6) = -4.289500e-01;
            coeffs(7) = -2.748200e-02;
            coeffs(8) = -4.168300e-01;
        }
        else if (j == 4) {
            coeffs(0) = -1.955400e-17;
            coeffs(1) = -7.373800e-17;
            coeffs(2) = 1.578900e-16;
            coeffs(3) = 3.917100e-01;
            coeffs(4) = 1.237500e-02;
            coeffs(5) = -4.039200e-01;
            coeffs(6) = 4.167300e-01;
            coeffs(7) = 1.316600e-02;
            coeffs(8) = -4.297200e-01;
        }
        else{
            cout << "ERROR! WRONG J IN GTO";
            exit(1);
        }

        basis.add_3_21G_ne(pos, coeffs); // for neon atom
        return basis.get_orb();

    }
    //else if (sys == "h") basis.add_3_21G_h(pos); // for hydrogen atom

    // delete from memory
    //basis.delContracted();
    //delete m_basis;

}

/*
void GTO::Contracted_orbital_1s()
{
    return basis.get_orb_1s();
}

void GTO::Contracted_orbital_2s()
{
    return basis.get_orb_2s();
}

void GTO::Contracted_orbital_3s()
{
    return basis.get_orb_3s();
}

void GTO::Contracted_orbital_1p()
{
    return basis.get_orb_1p();
}

void GTO::Contracted_orbital_2p()
{
    return basis.get_orb_2p();
}

*/
