#include "heliumsimplenumerically.h"
#include "trialfunction.h"
#include "vmcsolver.h"


#include <iostream>

using namespace std;

HeliumSimpleNumerically::HeliumSimpleNumerically()
{
//    cout << "Simple got created" << localEnergy << endl;
}

double HeliumSimpleNumerically::waveFunction(const mat &r, VMCSolver *solver)
{

    double argument = 0;
    for(int i = 0; i < solver->getNParticleshere(); i++) {
        double rSingleParticle = 0;
        for(int j = 0; j < solver->getNDimensions(); j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        argument += sqrt(rSingleParticle);
    }
    return exp(-argument * solver->getAlpha());
}
