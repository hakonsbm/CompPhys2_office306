#include "trialfunction.h"
#include "vmcsolver.h"

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;


TrialFunction::TrialFunction()
{

}

double TrialFunction::spinFactor(int i, int j)
{

    if(spin(i) == spin(j))
        return 1./4.;
    else
        return 1./2.;

}
