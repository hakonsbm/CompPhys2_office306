#include "vmcsolver.h"
#include <trialfunction.h>
#include "trialFunctions/heliumjastrowanalytical.h"
#include "trialFunctions/heliumjastrownumerical.h"
#include "trialFunctions/heliumsimpleanalytical.h"
#include "trialFunctions/heliumsimplenumerical.h"
#include "trialFunctions/hydrogen.h"
#include "trialFunctions/beryllium.h"
#include "trialFunctions/neon.h"
#include "trialFunctions/helium.h"
#include "trialFunctions/hydrogentwo.h"
#include "lib.h"

#include <unittest++/UnitTest++.h>
#include <mpi.h>
#include <armadillo>

using namespace arma;


TEST(Hydrogenic) {


//    VMCSolver *solver = new VMCSolver();

//    cout << endl << "Running Hydrogen test" << endl << endl;

//    solver->setTrialFunction(new Hydrogen(solver));
//    solver->runMonteCarloIntegration();
//    CHECK_EQUAL(0., solver->getEnergyVar());
//    CHECK_EQUAL(-1./2, solver->getEnergy());


//    cout << endl << "Running Helium test" << endl << endl;
//    solver->setTrialFunction(new Helium(solver));
//    solver->switchElectronInteraction(false);
//    solver->trialFunction()->setAnalytical(true);
//    solver->setAlpha(solver->getCharge());
//    solver->switchbBlockSampling(false);
//    solver->setCycles(10000);
//    solver->runMasterIntegration();
//    CHECK_EQUAL(0., solver->getEnergyVar());
//    CHECK_EQUAL(-4, solver->getEnergy());



//    cout << endl << "Running Beryllium test" << endl << endl;
//    solver->setTrialFunction(new Beryllium(solver));
//    solver->switchElectronInteraction(false);
//    solver->trialFunction()->setAnalytical(true);
//    solver->setAlpha(solver->getCharge());
//    solver->setCycles(10000);
//    solver->runMasterIntegration();
//    CHECK_EQUAL(0., solver->getEnergyVar());
//    CHECK_EQUAL(-20, solver->getEnergy());


//    cout << endl << "Running Neon test" << endl << endl;
//    solver->setTrialFunction(new Neon(solver));
//    solver->switchElectronInteraction(false);
//    solver->trialFunction()->setAnalytical(true);
//    solver->setCycles(1000);
//    solver->setAlpha(solver->getCharge());
//    solver->runMasterIntegration();
//    CHECK_EQUAL(0., solver->getEnergyVar());
//    CHECK_EQUAL(-200., solver->getEnergy());


}

TEST(HeliumJastrow)
{
    cout << "Testing The analytical version of the Helium wavefunction" << endl;

    double oldVersion, newVersion;


    //Testing the machinery for the gradient ratio of Psi_C vs a calcualted ratio
    VMCSolver *solver = new VMCSolver();

    solver->setTrialFunction(new Helium(solver));
//    solver->switchElectronInteraction(false);
//    solver->trialFunction()->setAnalytical(true);
    solver->setAlpha(solver->getCharge());
    solver->setBeta(2);
//    solver->setCycles(1000);

    mat r = zeros (2,3);
    vec correct = zeros (3);
    vec calculated = zeros(3);
    double r12 = 0;

    double particles = 2;
    double dimensions = 3;
    double beta = solver->getBeta();
    long idum = -1;

    for(int i = 0; i < particles; i ++ )
    {
        for(int j = 0; j < dimensions; j++)
        {
           r(i,j) = ran2(&idum);
        }
    }

    r12 = norm(r.row(0) - r.row(1));
    correct = ((r.row(0) - r.row(1)).t())  / (r12 * pow(1+(beta*r12), 2)) ;

    cout << solver->derivatives()->analyticalCorrelationGradient(r,solver)<< endl;



    CHECK_CLOSE( correct(0) , calculated(0), 0.001  );
    CHECK_CLOSE( correct(1) , calculated(1), 0.001  );
    CHECK_CLOSE( correct(2) , calculated(2), 0.001  );


    //Checking the laplacian ratio of Psi_C vs a calculated one.
//    cout << solver->derivatives()->analyticalCorrelationDoubleDerivative(r, solver) << endl;
    solver->setTrialFunction(new Helium(solver));
    solver->trialFunction()->setAnalytical(true);

}

TEST(FirstOrderDerivative)
{
//    VMCSolver *solver = new VMCSolver();

//    solver->setTrialFunction(new Helium(solver));


//    double particles = solver->getNParticles();
//    double dimensions = solver->getNDimensions();
//    long idum = -1;

//    mat r = zeros (particles,dimensions);
//    mat gradientNumerical = zeros(particles, dimensions);
//    mat gradientAnalytical = zeros(particles, dimensions);


//    //Random positions to test derivative
//    for(int i = 0; i < particles; i ++ )
//    {
//     for(int j = 0; j < dimensions; j++)
//     {
//        r(i,j) = ran2(&idum);
//     }
//    }
//    solver->trialFunction()->setAnalytical(false);
//    solver->determinant()->updateSlaterMatrices(r,solver);
//    solver->derivatives()->numericalGradient(gradientNumerical, r, solver);
////    //Testing
////    gradientNumerical = solver->determinant()->gradientSlaterDeterminant(r,solver);

//    solver->trialFunction()->setAnalytical(true);
//    solver->determinant()->updateSlaterMatrices(r,solver);
//    solver->derivatives()->analyticalGradient(gradientAnalytical, r , solver);
////    //Testing
////    gradientNumerical = solver->determinant()->gradientSlaterDeterminant(r,solver);

////    gradientAnalytical = solver->derivatives()->analyticalCorrelationGradient(r,solver);

//    if(solver->getRank() == 0)
//    {
//        cout << "End results" << endl;
//        cout << gradientNumerical << endl;
//        cout << gradientAnalytical << endl;
//    }

//    exit(0);

}

TEST(HYDROGENTWO_VS_HELIUM)
{
    //Make a test with R = 0 for hydrogenTwo where it should be the same as the Helium atom

//    cout << endl << "Running He vs H_2 test" << endl << endl;
//    VMCSolver *solver = new VMCSolver;
//    solver->setTrialFunction(new Helium(solver));
//    solver->trialFunction()->setAnalytical(false);
//    solver->switchbBlockSampling(false);
//    solver->setCycles(100000);
//    solver->runMasterIntegration();

//    mat r = zeros (2,3);
//    vec correct = zeros (3);
//    vec calculated = zeros(3);
//    double r12 = 0;

//    double particles = 2;
//    double dimensions = 3;
//    double beta = solver->getBeta();
//    long idum = -10000;

//    for(int i = 0; i < particles; i ++ )
//    {
//        for(int j = 0; j < dimensions; j++)
//        {
//           r(i,j) = ran2(&idum);
//        }
//    }

////    solver->trialFunction()->waveFunction(r, solver);
////    solver->trialFunction()->localEnergy(r,solver);
//    solver->determinant()->updateSlaterMatrices(r, solver);
//    cout << solver->trialFunction()->waveFunction(r,solver) << endl;


//    solver->setTrialFunction(new HydrogenTwo(solver));
//    solver->trialFunction()->setAnalytical(false);
//    solver->trialFunction()->setNucleusDistance(1.4);
////    solver->switchElectronInteraction(true);
//    solver->switchbBlockSampling(false);
//    solver->setCycles(100000);
//    solver->runMasterIntegration();



////    cout << r << endl;
//    solver->determinant()->updateSlaterMatrices(r, solver);
////    solver->trialFunction()->waveFunction(r, solver);
////    solver->trialFunction()->localEnergy(r,solver);
//    cout << solver->trialFunction()->waveFunction(r,solver) << endl;





}

TEST(AnalyticalHelium)
{


    //Testing the laplacian ratio


//    cout << "Testing The analytical version of the Helium wavefunction" << endl;

    double analytical, numerical, analyticalVar , numericalVar;
    VMCSolver *solver = new VMCSolver();


//    cout << endl << "Running Helium test" << endl << endl;
    solver->setTrialFunction(new Helium(solver));
    solver->switchElectronInteraction(true);
    solver->trialFunction()->setAnalytical(false);
    solver->switchbBlockSampling(false);
    solver->setCycles(1000000);
    solver->runMasterIntegration();

////TestStuff
//    double particles = solver->getNParticles();
//    double dimensions = solver->getNDimensions();
//    long idum = -clock();

//    mat r = zeros (particles,dimensions);

//    //Random positions to test derivative
//    for(int i = 0; i < particles; i ++ )
//    {
//     for(int j = 0; j < dimensions; j++)
//     {
//        r(i,j) = ran2(&idum);
//     }
//    }
//    numerical = solver->trialFunction()->localEnergy(r,solver);


    numerical = solver->getEnergy();
    numericalVar = solver->getEnergyVar();


    solver->setTrialFunction(new Helium(solver));
    solver->switchElectronInteraction(true);
    solver->trialFunction()->setAnalytical(true);
    solver->switchbBlockSampling(false);
    solver->setCycles(1000000);
    solver->runMasterIntegration();

//    analytical = solver->trialFunction()->localEnergy(r,solver);


    analytical = solver->getEnergy();
    analyticalVar = solver->getEnergyVar();


    if(solver->getRank() == 0)
    {
        cout << "Energy with the old one: " << numerical << " and " << numericalVar << endl;
        cout << "Energy with the new one: " << analytical << " and " << analyticalVar << endl;
    }


}

