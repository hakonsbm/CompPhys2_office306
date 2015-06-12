//This file is maintained by an external python script and should not be edited manually.
#include <armadillo>
#include "basisbank.h"
#include "contracted.h"

using namespace std;
using namespace arma;

basisbank::basisbank()
{
    //initContracted(new Contracted);
    Contracted contracted;
}



void basisbank::add_3_21G_be(const vec corePos, const vec c){
    // s-orbital
    contracted.add_primitive(71.88760000,0.06442630,0,0,0,corePos);
    contracted.add_primitive(10.72890000,0.36609600,0,0,0,corePos);
    contracted.add_primitive(2.22205000,0.69593400,0,0,0,corePos);
    contracted.contract_orb(c[0]);
    // s-orbital
    contracted.add_primitive(1.29548000,-0.42106400,0,0,0,corePos);
    contracted.add_primitive(0.26888100,1.22407000,0,0,0,corePos);
    contracted.contract_orb(c[1]);
    // s-orbital
    contracted.add_primitive(0.07735000,1.00000000,0,0,0,corePos);
    contracted.contract_orb(c[2]);
    // p-orbital
    contracted.add_primitive(1.29548000,0.20513200,1,0,0,corePos);
    contracted.add_primitive(0.26888100,0.88252800,1,0,0,corePos);
    contracted.contract_orb(c[3]);
    // p-orbital
    contracted.add_primitive(1.29548000,0.20513200,0,1,0,corePos);
    contracted.add_primitive(0.26888100,0.88252800,0,1,0,corePos);
    contracted.contract_orb(c[5]);
    // p-orbital
    contracted.add_primitive(1.29548000,0.20513200,0,0,1,corePos);
    contracted.add_primitive(0.26888100,0.88252800,0,0,1,corePos);
    contracted.contract_orb(c[7]);
    // p-orbital
    contracted.add_primitive(0.07735000,1.00000000,1,0,0,corePos);
    contracted.contract_orb(c[4]);
    // p-orbital
    contracted.add_primitive(0.07735000,1.00000000,0,1,0,corePos);
    contracted.contract_orb(c[6]);
    // p-orbital
    contracted.add_primitive(0.07735000,1.00000000,0,0,1,corePos);
    contracted.contract_orb(c[8]);
}



void basisbank::add_3_21G_ne(const vec corePos, const vec c){
    // s-orbital
    contracted.add_primitive(515.72400000,0.05814300,0,0,0,corePos);
    contracted.add_primitive(77.65380000,0.34795100,0,0,0,corePos);
    contracted.add_primitive(16.81360000,0.71071400,0,0,0,corePos);
    contracted.contract_orb(c[0]);
    // s-orbital
    contracted.add_primitive(12.48300000,-0.40992200,0,0,0,corePos);
    contracted.add_primitive(2.66451000,1.22431000,0,0,0,corePos);
    contracted.contract_orb(c[1]);
    // s-orbital
    contracted.add_primitive(0.60625000,1.00000000,0,0,0,corePos);
    contracted.contract_orb(c[2]);
    // p-orbital
    contracted.add_primitive(12.48300000,0.24746000,1,0,0,corePos);
    contracted.add_primitive(2.66451000,0.85174300,1,0,0,corePos);
    contracted.contract_orb(c[3]);
    // p-orbital
    contracted.add_primitive(12.48300000,0.24746000,0,1,0,corePos);
    contracted.add_primitive(2.66451000,0.85174300,0,1,0,corePos);
    contracted.contract_orb(c[5]);
    // p-orbital
    contracted.add_primitive(12.48300000,0.24746000,0,0,1,corePos);
    contracted.add_primitive(2.66451000,0.85174300,0,0,1,corePos);
    contracted.contract_orb(c[7]);
    // p-orbital
    contracted.add_primitive(0.60625000,1.00000000,1,0,0,corePos);
    contracted.contract_orb(c[4]);
    // p-orbital
    contracted.add_primitive(0.60625000,1.00000000,0,1,0,corePos);
    contracted.contract_orb(c[6]);
    // p-orbital
    contracted.add_primitive(0.60625000,1.00000000,0,0,1,corePos);
    contracted.contract_orb(c[8]);
}



void basisbank::add_3_21G_he(const vec corePos, const vec c){
    // s-orbital
    contracted.add_primitive(13.62670000,0.17523000,0,0,0,corePos);
    contracted.add_primitive(1.99935000,0.89348300,0,0,0,corePos);
    contracted.contract_orb(c[0]);
    // s-orbital
    contracted.add_primitive(0.38299300,1.00000000,0,0,0,corePos);
    contracted.contract_orb(c[1]);
}


double basisbank::get_orb() {return contracted.get_orb();}
