#include "heliumjastrowanalytical.h"


HeliumJastrowAnalytical::HeliumJastrowAnalytical(VMCSolver *solver)
{
    simpleFlag = false;
    m_outfileName = "HeliumJastrowAnalytical";

    solver->setCharge(2);
    solver->setNParticles(2);
    solver->setAlpha(1.843);
    solver->setBeta(0.34);


}

HeliumJastrowAnalytical::~HeliumJastrowAnalytical()
{

}

double HeliumJastrowAnalytical::waveFunction(const mat &r, VMCSolver *solver)
{
    double nParticles = solver->getNParticles();
    double nDimensions = solver->getNDimensions();
    double alpha = solver->getAlpha();
    double beta = solver->getBeta();


    //double r12;
    vec rpos(nParticles);
    for(int i = 0; i < nParticles; i++) {
        double rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        rpos[i] = sqrt(rSingleParticle);
    }
    // assuming 2 particles
    //   (ta fra elektron-elektron pot.)
    //r12 = abs(rpos[0] - rpos[1]);
    double r12 = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = i + 1; j < nParticles; j++) {
            r12 = 0;
            for(int k = 0; k < nDimensions; k++) {
                r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));

            }
            r12 = sqrt(r12);
        }
    }
//    cout << exp(r12 / (2.0*(1 + beta * r12))) << endl;
    return exp(-accu(rpos) * alpha) * exp(r12 / (2.0*(1 + beta * r12))) ;
}

double HeliumJastrowAnalytical::localEnergy(const mat &r, VMCSolver *solver)
{
    double localEnergy1 = 0;
    double localEnergy2 = 0;

    //Grabbing all the necessary constants stored in the solver
    double charge = solver->getCharge();
    double alpha = solver->getAlpha();
    double beta = solver->getBeta();


    //Precalculating the distances
    double r1 = norm(r.row(0));
    double r2 = norm(r.row(1));
    double r12 = norm(r.row(0) - r.row(1));

    //The part of the local Energy that was the same as for the first simple case
    localEnergy1 = (alpha - charge)*(1./r1+1./r2) + 1./(r12) - pow(alpha,2);

    //The second contribution to the local Energy, I think there is an error here
    localEnergy2 = 0.5*pow((1+beta*r12),-2) * (alpha*(r1+r2)/r12 * (1. - norm_dot(r.row(0),r.row(1))) -
                                                    0.5*pow((1+beta*r12),-2) - 2./r12 + 2.*beta/(1 + beta * r12) );
//    cout <<  endl << "In JastrowAnalytical" << endl;
//    cout << "LocalEnergy1 is " <<localEnergy1 << endl;
//    cout << "LocalEnergy2 is " <<localEnergy2 << endl;

    return localEnergy1 + localEnergy2;
}

double HeliumJastrowAnalytical::lnDerivativeWaveFunction(const mat &r, VMCSolver *solver)
{
    return 0;
}

double HeliumJastrowAnalytical::lnSecondDerivativeWaveFunction(const mat &r, VMCSolver *solver)
{
    return 0;
}

