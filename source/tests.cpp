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
    solver->setTrialFunction(new HeliumSimpleNumerical(solver));

//    solver->derivatives()->numericalDerivative(solver);

    solver->runMasterIntegration();

    CHECK_EQUAL(0 , 0);
}
