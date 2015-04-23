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
    double alpha = solver->getAlpha();
    vec rpos(solver->getNParticles());
    double SD;


    for(int i = 0; i < solver->getNParticles(); i++) {
        double rSingleParticle = 0;
        for(int j = 0; j < solver->getNDimensions(); j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        rpos[i] = sqrt(rSingleParticle);
    }

    SD = solver->determinant()->calculateDeterminant(r,alpha,solver);


//    cout << "Starting new run" << endl << endl;
//    cout << "detUp is" << solver->determinant()->detUp << endl;
//    cout << "detDown is" << solver->determinant()->detDown << endl;

//    cout << "SlaterDeterminant is " << SD << endl;
//    cout << "Correct wavefunction is: "  << exp(-accu(rpos) * alpha) << endl << endl;

    return SD;//exp(-accu(rpos) * alpha);
}

double HeliumSimpleAnalytical::localEnergy(const mat &r, VMCSolver *solver)
{
    double kineticEnergy, potentialEnergy;

    double r1 = norm(r.row(0));
    double r2 = norm(r.row(1));
    double r12 = norm(r.row(0) - r.row(1));

    double charge = solver->getCharge();

//    solver->determinant()->updateSlaterMatrices(r,solver);

    kineticEnergy = solver->determinant()->laplacianSlaterDeterminant(r,solver)/(-2.);

//    kineticEnergy = solver->derivatives()->analyticalSimpleDoubleDerivative(r,solver)/(-2.);
//    cout << "r: " << endl << r << endl;
//    cout << "r1: " << r1 << endl;
//    cout << "r2: " << r2 << endl;

//    cout << "Attempt at new: " << solver->determinant()->laplacianSlaterDeterminant(r,solver)/(-2.) << endl;
//    cout << "Correct Energy: " << kineticEnergy << endl;


//    kineticEnergy = alpha*(1./r1+1./r2) - pow(alpha,2);

    //Taking away the electron electron interaction, used for some tests with Hydrogenic wavesfunctions
    if(solver->getElectronInteration())
        potentialEnergy = - charge*(1./r1+1./r2) + 1./(r12);
    else
        potentialEnergy = - charge*(1./r1+1./r2);


    //Returns the local energy, EL = (a-Z)(1/r1+1/r2)+1/r12-alpha^2)
        return kineticEnergy + potentialEnergy;

}
