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
    void add_3_21G_be(const vec corePos, const vec c);
    void add_3_21G_ne(const vec corePos, const vec c);
    void add_3_21G_he(const vec corePos, const vec c);
    //void initContracted(Contracted *contracted) {m_contracted = contracted;}
    //void delContracted() {delete m_contracted;}
    //Contracted *contracted() {return m_contracted;}
    double get_orb();
private:
    Contracted contracted;
};
#endif // BASISBANK_H
