#include "trialfunction.h"
#include "vmcsolver.h"

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;


TrialFunction::TrialFunction()
{
    m_nucleusDistance = 0;
    m_conjugateMethod = false ;
}

double TrialFunction::spinFactor(int i, int j)
{

    if(spin(i) == spin(j))
        return 1./4.;
    else
        return 1./2.;

}


