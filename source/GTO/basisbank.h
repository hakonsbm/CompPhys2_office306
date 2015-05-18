//This file is maintained by an external python script and should not be edited manually.
#ifndef BASISBANK_H
#define BASISBANK_H
#include <armadillo>
#include <contracted.h>

using namespace std;
using namespace arma;
 
class basisbank{
public:
    
    void add_3_21G_ne(const vec corePos);
    void add_3_21G_he(const vec corePos);
    void add_3_21G_h(const vec corePos);
    double get_contracted() {return Contracted->get_contracted();}
};
#endif // BASISBANK_H