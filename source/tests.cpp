#include "vmcsolver.h"
#include <trialfunction.h>
#include "trialFunctions/heliumjastrowanalytical.h"
#include "trialFunctions/heliumjastrownumerical.h"
#include "trialFunctions/heliumsimpleanalytical.h"
#include "trialFunctions/heliumsimplenumerical.h"
#include "trialFunctions/hydrogen.h"
#include "trialFunctions/beryllium.h"
#include "trialFunctions/neon.h"

#include <unittest++/UnitTest++.h>
#include <mpi.h>


TEST(Hydrogen) {

    cout << endl << "Running Hydrogen test" << endl << endl;

    VMCSolver *solver = new VMCSolver();
    solver->setTrialFunction(new Hydrogen(solver));
    solver->runMonteCarloIntegration();
    CHECK_EQUAL(0., solver->getEnergyVar());

}

TEST(Derivatives)
{
    cout << endl << "Running Derivatives test" << endl << endl;
    VMCSolver *solver = new VMCSolver();
    solver->setTrialFunction(new HeliumSimpleAnalytical(solver));

//    solver->derivatives()->numericalDerivative(solver);

    double analytic, numerical;

    solver->runMasterIntegration();

    analytic = solver->getEnergy();


    solver->setTrialFunction(new HeliumSimpleNumerical(solver));

    solver->runMasterIntegration();

    numerical = solver->getEnergy();

    cout << "The difference was " << analytic - numerical << endl;



    MPI_Finalize ();

    CHECK_CLOSE(fabs(analytic - numerical) , 0, 0.05 );
}
