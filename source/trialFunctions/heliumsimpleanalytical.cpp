#include "heliumsimpleanalytical.h"

HeliumSimpleAnalytical::HeliumSimpleAnalytical(VMCSolver *solver)
{
    simpleFlag = true;
    m_outfileName = "HeliumSimpleAnalytical";

    solver->setCharge(2);
    solver->setNParticles(2);
    solver->setAlpha(1.65);
    solver->setBeta(0);
}

double HeliumSimpleAnalytical::waveFunction(const mat &r, VMCSolver *solver)
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

double HeliumSimpleAnalytical::localEnergy(const mat &r, VMCSolver *solver)
{
    double kineticEnergy, potentialEnergy;

    double r1 = norm(r.row(0));
    double r2 = norm(r.row(1));
    double r12 = norm(r.row(0) - r.row(1));

    double alpha = solver->getAlpha();
    double charge = solver->getCharge();

    kineticEnergy = solver->derivatives()->analyticalSimpleDoubleDerivative(r,solver);
//    kineticEnergy = alpha*(1./r1+1./r2) - pow(alpha,2);

    //Taking away the electron electron interaction, used for some tests with Hydrogenic wavesfunctions
    if(solver->getElectronInteration())
        potentialEnergy = - charge*(1./r1+1./r2) + 1./(r12);
    else
        potentialEnergy = - charge*(1./r1+1./r2);


    //Returns the local energy, EL = (a-Z)(1/r1+1/r2)+1/r12-alpha^2)
        return kineticEnergy + potentialEnergy;

}
