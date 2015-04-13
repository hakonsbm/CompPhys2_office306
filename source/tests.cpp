#include "vmcsolver.h"
#include <trialfunction.h>
#include "trialFunctions/heliumjastrowanalytical.h"
#include "trialFunctions/heliumjastrownumerical.h"
#include "trialFunctions/heliumsimpleanalytical.h"
#include "trialFunctions/heliumsimplenumerical.h"
#include "trialFunctions/hydrogen.h"
#include "trialFunctions/beryllium.h"
#include "trialFunctions/neon.h"
#include "lib.h"

#include <unittest++/UnitTest++.h>
#include <mpi.h>
#include <armadillo>

using namespace arma;


TEST(Hydrogenic) {


    cout << endl << "Running Hydrogen test" << endl << endl;

    VMCSolver *solver = new VMCSolver();
    solver->setTrialFunction(new Hydrogen(solver));
    solver->runMonteCarloIntegration();
    CHECK_EQUAL(0., solver->getEnergyVar());
    CHECK_EQUAL(-1./2, solver->getEnergy());


    cout << endl << "Running Helium test" << endl << endl;
    solver->setTrialFunction(new HeliumSimpleAnalytical(solver));
    solver->switchElectronInteraction(false);
    solver->setAlpha(solver->getCharge());
    solver->runMasterIntegration();
    CHECK_EQUAL(0., solver->getEnergyVar());
    CHECK_EQUAL(-4, solver->getEnergy());
}


TEST(Derivatives)
{
    cout << endl << "Running Derivatives test" << endl << endl;


    VMCSolver *solver = new VMCSolver();
    solver->setTrialFunction(new HeliumSimpleAnalytical(solver));

//    double nParticles = solver->getNParticles();
//    double nDimensions = solver->getNDimensions();
//    long int idum = -1;
//    mat r = zeros<mat>(nParticles, nDimensions);

//    //random position to test the analytical derivation against the numerical;
//    for(int i = 0; i < nParticles; i++) {
//        for(int j = 0; j < nDimensions(); j++) {
//            r(i,j) = GaussianDeviate(&idum)*sqrt(stepLength);
//        }
//    }

//////    solver->derivatives()->numericalDerivative(solver);

    double analytic, numerical;
//    analytic = numerical;


    solver->runMasterIntegration();

    analytic = solver->getEnergy();


    solver->setTrialFunction(new HeliumSimpleNumerical(solver));

    solver->runMasterIntegration();

    numerical = solver->getEnergy();

    cout << "The difference was " << analytic - numerical << endl;




    CHECK_CLOSE(fabs(analytic - numerical) , 0, 0.05 );

    MPI_Finalize ();

}
