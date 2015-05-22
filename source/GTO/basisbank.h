//This file is maintained by an external python script and should not be edited manually.
#ifndef BASISBANK_H
#define BASISBANK_H
#include <armadillo>
#include "contracted.h"

using namespace std;
using namespace arma;
 
class basisbank{
public:
    basisbank();
    void add_3_21G_be(const vec corePos);
    void add_3_21G_ne(const vec corePos);
    void add_3_21G_he(const vec corePos);
    void add_3_21G_h(const vec corePos);
    void initContracted(Contracted *contracted) {m_contracted = contracted;}
    Contracted *contracted() {return m_contracted;}
    double get_orb_1s();
    double get_orb_2s();
    double get_orb_pX();
    double get_orb_pY();
    double get_orb_pZ();
private:
    Contracted *m_contracted;
};
#endif // BASISBANK_H