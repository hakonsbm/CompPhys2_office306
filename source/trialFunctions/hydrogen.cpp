#include "hydrogen.h"

Hydrogen::Hydrogen(VMCSolver* solver)
{
    simpleFlag = true;
    m_outfileName = "Hydrogen";

    solver->setCharge(1);
    solver->setNParticles(1);
    solver->setAlpha(1.);
    solver->setBeta(0);
}

Hydrogen::~Hydrogen()
{

}

double Hydrogen::waveFunction(const mat &r, VMCSolver *solver)
{
    vec rpos(solver->getNParticles());
    for(int i = 0; i < solver->getNParticles(); i++) {
        double rSingleParticle = 0;
        for(int j = 0; j < solver->getNDimensions(); j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        rpos[i] = sqrt(rSingleParticle);
    }
    return exp(-accu(rpos) * solver->getAlpha());
}

double Hydrogen::localEnergy(const mat &r, VMCSolver *solver)
{

    double r1 = norm(r.row(0));


    double alpha = solver->getAlpha();
    double charge = solver->getCharge();

    //Returns the local energy, EL = (a-Z)(1/r1+1/r2)+1/r12-alpha^2)
    return -pow(r1,-1) - (alpha/2)* (alpha - (2/r1)) ;
}
