//This file is maintained by an external python script and should not be edited manually.
#include <basisbank.h>
#include <armadillo>
#include <basis.h>
#include <primitive.h>
using namespace std;
using namespace arma;
basisbank::basisbank(basis BS){
    bs = BS;} 
basisbank::basisbank(){} 
 
//#  3-21G  EMSL  Basis Set Exchange Library  5/13/15 3:46 AM

//# Elements                             References

//# --------                             ----------

//#  H - Ne: J.S. Binkley, J.A. Pople, W.J. Hehre, J. Am. Chem. Soc 102 939 (1980)

//# Na - Ar: M.S. Gordon, J.S. Binkley, J.A. Pople, W.J. Pietro and W.J. Hehre, 

//#          J. Am. Chem. Soc. 104, 2797 (1983).

//#  K - Ca: K.D. Dobbs, W.J. Hehre, J. Comput. Chem. 7, 359 (1986). 

//# Ga - Kr: K.D. Dobbs, W.J. Hehre, J. Comput. Chem. 7, 359 (1986).

//# Sc - Zn: K.D. Dobbs, W.J. Hehre, J. Comput. Chem. 8, 861 (1987). 

//#  Y - Cd: K.D. Dobbs, W.J. Hehre, J. Comput. Chem. 8, 880 (1987). 

//# Cs     : A 3-21G quality set derived from the Huzinage MIDI basis sets.

//#          E.D. Glendening and D. Feller, J. Phys. Chem. 99, 3060 (1995)

//#   

void basisbank::add_3_21G_ne(vec3 corePos){

    bs.add_state();
    Primitive S0A0 = bs.turbomolePrimitive(0.05814300,515.72400000,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A0);
    Primitive S1A0 = bs.turbomolePrimitive(0.34795100,77.65380000,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A0);
    Primitive S2A0 = bs.turbomolePrimitive(0.71071400,16.81360000,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A0);
    bs.add_state();
    Primitive S0A1 = bs.turbomolePrimitive(-0.40992200,12.48300000,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A1);
    Primitive S1A1 = bs.turbomolePrimitive(1.22431000,2.66451000,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A1);
    bs.add_state();
    Primitive S0A2 = bs.turbomolePrimitive(1.00000000,0.60625000,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A2);
    bs.add_state();
    Primitive P0A3 = bs.turbomolePrimitive(0.24746000,12.48300000,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0A3);
    Primitive P0B3 = bs.turbomolePrimitive(0.24746000,12.48300000,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0B3);
    Primitive P0C3 = bs.turbomolePrimitive(0.24746000,12.48300000,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0C3);
    Primitive P1A3 = bs.turbomolePrimitive(0.85174300,2.66451000,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1A3);
    Primitive P1B3 = bs.turbomolePrimitive(0.85174300,2.66451000,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1B3);
    Primitive P1C3 = bs.turbomolePrimitive(0.85174300,2.66451000,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1C3);
    bs.add_state();
    Primitive P0A4 = bs.turbomolePrimitive(1.00000000,0.60625000,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0A4);
    Primitive P0B4 = bs.turbomolePrimitive(1.00000000,0.60625000,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0B4);
    Primitive P0C4 = bs.turbomolePrimitive(1.00000000,0.60625000,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0C4);
}


//#  3-21G  EMSL  Basis Set Exchange Library  5/13/15 4:09 AM

//# Elements                             References

//# --------                             ----------

//#  H - Ne: J.S. Binkley, J.A. Pople, W.J. Hehre, J. Am. Chem. Soc 102 939 (1980)

//# Na - Ar: M.S. Gordon, J.S. Binkley, J.A. Pople, W.J. Pietro and W.J. Hehre, 

//#          J. Am. Chem. Soc. 104, 2797 (1983).

//#  K - Ca: K.D. Dobbs, W.J. Hehre, J. Comput. Chem. 7, 359 (1986). 

//# Ga - Kr: K.D. Dobbs, W.J. Hehre, J. Comput. Chem. 7, 359 (1986).

//# Sc - Zn: K.D. Dobbs, W.J. Hehre, J. Comput. Chem. 8, 861 (1987). 

//#  Y - Cd: K.D. Dobbs, W.J. Hehre, J. Comput. Chem. 8, 880 (1987). 

//# Cs     : A 3-21G quality set derived from the Huzinage MIDI basis sets.

//#          E.D. Glendening and D. Feller, J. Phys. Chem. 99, 3060 (1995)

//#   

void basisbank::add_3_21G_he(vec3 corePos){

    bs.add_state();
    Primitive S0A0 = bs.turbomolePrimitive(0.17523000,13.62670000,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A0);
    Primitive S1A0 = bs.turbomolePrimitive(0.89348300,1.99935000,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A0);
    bs.add_state();
    Primitive S0A1 = bs.turbomolePrimitive(1.00000000,0.38299300,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A1);
}


//# 

void basisbank::add_3_21G_h(vec3 corePos){

    bs.add_state();
    Primitive S0A0 = bs.turbomolePrimitive(0.15628500,5.44717800,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A0);
    Primitive S1A0 = bs.turbomolePrimitive(0.90469100,0.82454700,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A0);
    bs.add_state();
    Primitive S0A1 = bs.turbomolePrimitive(1.00000000,0.18319200,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A1);
}


