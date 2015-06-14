#include "trialfunction.h"
#include "vmcsolver.h"

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;


TrialFunction::TrialFunction()
{
//    m_nucleusDistance = 1.4;
    m_nucleusDistance = 4.63;

    m_conjugateMethod = false ;
}

double TrialFunction::spinFactor(int i, int j)
{

    if(spin(i) == spin(j))
        return 1./4.;
    else
        return 1./2.;

}

void TrialFunction::setNucleusDistance(double R)
{

    m_nucleusDistance = R;
    if(R < 0.001)
    {
        m_zeroDistance = true;
    }
    else m_zeroDistance = false;
}

