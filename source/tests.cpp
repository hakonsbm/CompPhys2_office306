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


    VMCSolver *solver = new VMCSolver();

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


TEST(Gradients)
{
//    VMCSolver *solver = new VMCSolver();

//    solver->setTrialFunction(new Helium(solver));

//    double particles = solver->getNParticles();
//    double dimensions = solver->getNDimensions();
//    long idum = clock();

//    mat r = zeros (particles,dimensions);
//    mat gradientNumerical = zeros(particles, dimensions);
//    mat gradientAnalytical = zeros(particles, dimensions);


//    //Random positions to test gradient
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


//    solver->trialFunction()->setAnalytical(true);
//    solver->determinant()->updateSlaterMatrices(r,solver);
//    solver->derivatives()->analyticalGradient(gradientAnalytical, r , solver);


//    if(solver->getRank() == 0)
//    {
//        cout << "End results" << endl;
//        cout << gradientNumerical << endl;
//        cout << gradientAnalytical << endl;
//    }


}

TEST(AnalyticalHelium)
{
//    //This tests the analytical machinery calculating the local energy by testing it against
//    // Helium for which we have a easy closed solution. This tests that Helium using the machinery gets the same E_L value with a random
//    // electron position
//    double analytical, numerical;
//    VMCSolver *solver = new VMCSolver();


//    solver->setTrialFunction(new HeliumJastrowAnalytical(solver));
//    solver->switchElectronInteraction(true);
//    solver->trialFunction()->setAnalytical(true);


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

//    solver->determinant()->updateSlaterMatrices(r,solver);
//    numerical = solver->trialFunction()->localEnergy(r,solver);


//    solver->setTrialFunction(new Helium(solver));
//    solver->switchElectronInteraction(true);
//    solver->trialFunction()->setAnalytical(true);


//    solver->determinant()->updateSlaterMatrices(r,solver);
//    analytical = solver->trialFunction()->localEnergy(r,solver);

//    if(solver->getRank() == 0)
//    {
//        cout << "Energy with the correct E_L one: " << numerical << endl;
//        cout << "Energy with the E_L machinery : " << analytical << endl;
//    }


//    CHECK_CLOSE(numerical, analytical, 0.00001 );
}


TEST(HYDROGENTWO_VS_HELIUM)
{
//    //Make a test with R = 0 for hydrogenTwo where it should be the same as the Helium atom

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

TEST(ANALYTICAL_VS_NUMERICAL)
{
    double analytical, numerical;
    VMCSolver *solver = new VMCSolver();


//    solver->setTrialFunction(new Beryllium(solver));
//    solver->switchElectronInteraction(true);
//    solver->trialFunction()->setAnalytical(false);
//    solver->setCycles(50000);
//    solver->runMasterIntegration();

//    numerical = solver->getEnergy();

//    solver->setTrialFunction(new Beryllium(solver));
//    solver->switchElectronInteraction(true);
//    solver->trialFunction()->setAnalytical(true);
//    solver->setCycles(50000);
//    solver->runMasterIntegration();

//    analytical = solver->getEnergy();

//    if(solver->getRank() == 0)
//    {
//        cout << "Energy with the correct E_L one: " << numerical << endl;
//        cout << "Energy with the E_L machinery : " << analytical << endl;
//    }


//    CHECK_CLOSE(numerical, analytical, 0.01 );
//    for(int i = 0; i < 100500; i++)
    {
        cout << "Running Neon numerical" << endl;
        solver->setTrialFunction(new Neon(solver));
        solver->switchElectronInteraction(true);
        solver->trialFunction()->setAnalytical(false);
        solver->setCycles(10000);
//        solver->runMasterIntegration();

        numerical = solver->getEnergy();

//    //TestStuff
//        double particles = solver->getNParticles();
//        double dimensions = solver->getNDimensions();
//        long idum = -clock();

//        mat r = zeros (particles,dimensions);

//        //Random positions to test derivative
//        for(int i = 0; i < particles; i ++ )
//        {
//         for(int j = 0; j < dimensions; j++)
//         {
//            r(i,j) = ran2(&idum);

//         }
//        }

//        solver->determinant()->updateSlaterMatrices(r,solver);
//        numerical = solver->trialFunction()->localEnergy(r,solver);


        cout << "Running Neon analytical" << endl;

        solver->setTrialFunction(new Beryllium(solver));
        solver->switchElectronInteraction(true);
        solver->trialFunction()->setAnalytical(true);
        solver->setCycles(100000);
        solver->runMasterIntegration();

        analytical = solver->getEnergy();

//        solver->determinant()->updateSlaterMatrices(r,solver);
//        analytical = solver->trialFunction()->localEnergy(r,solver);


        if(solver->getRank() == 0)
        {
            cout << "Energy with the correct E_L one: " << numerical << endl;
            cout << "Energy with the E_L machinery : " << analytical << endl;
        }

//        cout << solver->determinant()->detUpInverseOld <<endl;

//        if(fabs(numerical - analytical) > 10000)
//            break;
    }


    CHECK_CLOSE(numerical, analytical, 0.01 );
}


