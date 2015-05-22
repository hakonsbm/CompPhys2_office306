//This file is maintained by an external python script and should not be edited manually.
#include <armadillo>
#include "basisbank.h"
#include "contracted.h"

using namespace std;
using namespace arma;

basisbank::basisbank()
{
    initContracted(new Contracted);
}



void basisbank::add_3_21G_be(const vec corePos){
    // s-orbital
    contracted()->add_primitive(71.88760000,0.06442630,0,0,0,corePos);
    contracted()->add_primitive(10.72890000,0.36609600,0,0,0,corePos);
    contracted()->add_primitive(2.22205000,0.69593400,0,0,0,corePos);
    contracted()->contract_orb_2s();
    // s-orbital
    contracted()->add_primitive(1.29548000,-0.42106400,0,0,0,corePos);
    contracted()->add_primitive(0.26888100,1.22407000,0,0,0,corePos);
    contracted()->contract_orbK();
    // s-orbital
    contracted()->add_primitive(0.07735000,1.00000000,0,0,0,corePos);
    contracted()->contract_orb_1s();
    // p-orbital
    contracted()->add_primitive(1.29548000,0.20513200,1,0,0,corePos);
    contracted()->add_primitive(0.26888100,0.88252800,1,0,0,corePos);
    contracted()->contract_orb_pX();
    // p-orbital
    contracted()->add_primitive(1.29548000,0.20513200,0,1,0,corePos);
    contracted()->add_primitive(0.26888100,0.88252800,0,1,0,corePos);
    contracted()->contract_orb_pY();
    // p-orbital
    contracted()->add_primitive(1.29548000,0.20513200,0,0,1,corePos);
    contracted()->add_primitive(0.26888100,0.88252800,0,0,1,corePos);
    contracted()->contract_orb_pZ();
    // p-orbital
    contracted()->add_primitive(0.07735000,1.00000000,1,0,0,corePos);
    contracted()->contract_orb_pX();
    // p-orbital
    contracted()->add_primitive(0.07735000,1.00000000,0,1,0,corePos);
    contracted()->contract_orb_pY();
    // p-orbital
    contracted()->add_primitive(0.07735000,1.00000000,0,0,1,corePos);
    contracted()->contract_orb_pZ();
}



void basisbank::add_3_21G_ne(const vec corePos){
    // s-orbital
    contracted()->add_primitive(515.72400000,0.05814300,0,0,0,corePos);
    contracted()->add_primitive(77.65380000,0.34795100,0,0,0,corePos);
    contracted()->add_primitive(16.81360000,0.71071400,0,0,0,corePos);
    contracted()->contract_orb_2s();
    // s-orbital
    contracted()->add_primitive(12.48300000,-0.40992200,0,0,0,corePos);
    contracted()->add_primitive(2.66451000,1.22431000,0,0,0,corePos);
    contracted()->contract_orbK();
    // s-orbital
    contracted()->add_primitive(0.60625000,1.00000000,0,0,0,corePos);
    contracted()->contract_orb_1s();
    // p-orbital
    contracted()->add_primitive(12.48300000,0.24746000,1,0,0,corePos);
    contracted()->add_primitive(2.66451000,0.85174300,1,0,0,corePos);
    contracted()->contract_orb_pX();
    // p-orbital
    contracted()->add_primitive(12.48300000,0.24746000,0,1,0,corePos);
    contracted()->add_primitive(2.66451000,0.85174300,0,1,0,corePos);
    contracted()->contract_orb_pY();
    // p-orbital
    contracted()->add_primitive(12.48300000,0.24746000,0,0,1,corePos);
    contracted()->add_primitive(2.66451000,0.85174300,0,0,1,corePos);
    contracted()->contract_orb_pZ();
    // p-orbital
    contracted()->add_primitive(0.60625000,1.00000000,1,0,0,corePos);
    contracted()->contract_orb_pX();
    // p-orbital
    contracted()->add_primitive(0.60625000,1.00000000,0,1,0,corePos);
    contracted()->contract_orb_pY();
    // p-orbital
    contracted()->add_primitive(0.60625000,1.00000000,0,0,1,corePos);
    contracted()->contract_orb_pZ();
}



void basisbank::add_3_21G_he(const vec corePos){
    // s-orbital
    contracted()->add_primitive(13.62670000,0.17523000,0,0,0,corePos);
    contracted()->add_primitive(1.99935000,0.89348300,0,0,0,corePos);
    contracted()->contract_orbK();
    // s-orbital
    contracted()->add_primitive(0.38299300,1.00000000,0,0,0,corePos);
    contracted()->contract_orb_1s();
}



void basisbank::add_3_21G_h(const vec corePos){
    // s-orbital
    contracted()->add_primitive(5.44717800,0.15628500,0,0,0,corePos);
    contracted()->add_primitive(0.82454700,0.90469100,0,0,0,corePos);
    contracted()->contract_orbK();
    // s-orbital
    contracted()->add_primitive(0.18319200,1.00000000,0,0,0,corePos);
    contracted()->contract_orb_1s();
}


double basisbank::get_orb_1s() {return contracted()->get_orb_1s();}
double basisbank::get_orb_2s() {return contracted()->get_orb_2s();}
double basisbank::get_orb_pX() {return contracted()->get_orb_pX();}
double basisbank::get_orb_pY() {return contracted()->get_orb_pY();}
double basisbank::get_orb_pZ() {return contracted()->get_orb_pZ();}
