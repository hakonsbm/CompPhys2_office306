//This file is maintained by an external python script and should not be edited manually.
#include <armadillo>
#include <basisbank.h>
#include <contracted.h>

using namespace std;
using namespace arma;

basisbank::basisbank()
{
contracted Contracted = new contracted();
} 


void basisbank::add_3_21G_ne(const vec corePos){
    // s-orbital
    Contracted add_primitive(0.05814300,515.72400000,0,0,0,corePos);
    Contracted add_primitive(0.34795100,77.65380000,0,0,0,corePos);
    Contracted add_primitive(0.71071400,16.81360000,0,0,0,corePos);
    // s-orbital
    Contracted add_primitive(-0.40992200,12.48300000,0,0,0,corePos);
    Contracted add_primitive(1.22431000,2.66451000,0,0,0,corePos);
    // s-orbital
    Contracted add_primitive(1.00000000,0.60625000,0,0,0,corePos);
    // p-orbital
    Contracted add_primitive(0.24746000,12.48300000,1,0,0,corePos);
    Contracted add_primitive(0.24746000,12.48300000,0,1,0,corePos);
    Contracted add_primitive(0.24746000,12.48300000,0,0,1,corePos);
    Contracted add_primitive(0.85174300,2.66451000,1,0,0,corePos);
    Contracted add_primitive(0.85174300,2.66451000,0,1,0,corePos);
    Contracted add_primitive(0.85174300,2.66451000,0,0,1,corePos);
    // p-orbital
    Contracted add_primitive(1.00000000,0.60625000,1,0,0,corePos);
    Contracted add_primitive(1.00000000,0.60625000,0,1,0,corePos);
    Contracted add_primitive(1.00000000,0.60625000,0,0,1,corePos);
}



void basisbank::add_3_21G_he(const vec corePos){
    // s-orbital
    Contracted add_primitive(0.17523000,13.62670000,0,0,0,corePos);
    Contracted add_primitive(0.89348300,1.99935000,0,0,0,corePos);
    // s-orbital
    Contracted add_primitive(1.00000000,0.38299300,0,0,0,corePos);
}



void basisbank::add_3_21G_h(const vec corePos){
    // s-orbital
    Contracted add_primitive(0.15628500,5.44717800,0,0,0,corePos);
    Contracted add_primitive(0.90469100,0.82454700,0,0,0,corePos);
    // s-orbital
    Contracted add_primitive(1.00000000,0.18319200,0,0,0,corePos);
}


