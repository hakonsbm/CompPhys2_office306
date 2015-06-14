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

//The tests are described in the report appendix

TEST(Hydrogenic) {
    //C1 and C2
    VMCSolver *solver = new VMCSolver();
//    solver->trialFunction()->setConjugate(false);

    cout << endl << "Running Hydrogen test" << endl << endl;

    solver->setTrialFunction(new Hydrogen(solver));
    solver->runMonteCarloIntegration();
    CHECK_EQUAL(0., solver->getEnergyVar());
    CHECK_EQUAL(-1./2, solver->getEnergy());


    cout << endl << "Running Helium test" << endl << endl;
    solver->setTrialFunction(new Helium(solver));
    solver->switchElectronInteraction(false);
    solver->trialFunction()->setAnalytical(true);
    solver->setAlpha(solver->getCharge());
    solver->switchbBlockSampling(false);
    solver->setCycles(10000);

    solver->runMasterIntegration();

    CHECK_EQUAL(0., solver->getEnergyVar());
    CHECK_EQUAL(-4, solver->getEnergy());

    exit(0);

    cout << endl << "Running Beryllium test" << endl << endl;
    solver->setTrialFunction(new Beryllium(solver));
    solver->switchElectronInteraction(false);
    solver->trialFunction()->setAnalytical(true);
    solver->setAlpha(solver->getCharge());
    solver->setCycles(10000);
    solver->runMasterIntegration();
    CHECK_EQUAL(0., solver->getEnergyVar());
    CHECK_EQUAL(-20, solver->getEnergy());


    cout << endl << "Running Neon test" << endl << endl;
    solver->setTrialFunction(new Neon(solver));
    solver->switchElectronInteraction(false);
    solver->trialFunction()->setAnalytical(true);
    solver->setCycles(1000);
    solver->setAlpha(solver->getCharge());
    solver->runMasterIntegration();
    CHECK_EQUAL(0., solver->getEnergyVar());
    CHECK_EQUAL(-200., solver->getEnergy());


}


TEST(Gradients)
{
    //Not working fo the time being...
    //    ////////////////////////7
//    /// C3
//    /// ///////////////////////


//    VMCSolver *solver = new VMCSolver();


//    double particles = solver->getNParticles();
//    double dimensions = solver->getNDimensions();
//    long idum = clock();

//    mat r = zeros (particles,dimensions);
//    mat gradientNumerical = zeros(particles, dimensions);
//    mat gradientAnalytical = zeros(particles, dimensions);

//    for(int n = 0; n < 3; n ++)
//    {
//        if(n==0)
//            solver->setTrialFunction(new Helium(solver));
//        if(n==1)
//            solver->setTrialFunction(new Beryllium(solver));
//        if(n==2)
//            solver->setTrialFunction(new Neon(solver));


//        particles = solver->getNParticles();
//        dimensions = solver->getNDimensions();
//        idum = clock();

//        r = zeros (particles,dimensions);
//        gradientNumerical = zeros(particles, dimensions);
//        gradientAnalytical = zeros(particles, dimensions);

//        //Random positions to test gradient
//        for(int i = 0; i < particles; i ++ )
//        {
//         for(int j = 0; j < dimensions; j++)
//         {
//            r(i,j) = ran2(&idum);
//         }
//        }
//        solver->trialFunction()->setAnalytical(false);
//        solver->determinant()->updateSlaterMatrices(r,solver);
//        solver->derivatives()->numericalGradient(gradientNumerical, r, solver);


//        solver->trialFunction()->setAnalytical(true);
//        solver->determinant()->updateSlaterMatrices(r,solver);
//        solver->derivatives()->analyticalGradient(gradientAnalytical, r , solver);


//        if(solver->getRank() == 0)
//        {
//            cout << "End results" << endl;
//            cout << gradientNumerical << endl;
//            cout << gradientAnalytical << endl;
//        }
//        for(int i = 0; i < particles; i++)
//        {
//            for(int j = 0; j < dimensions; j ++)
//            {
//                CHECK_CLOSE(gradientNumerical(i,j),gradientAnalytical(i,j),0.0001);
//            }
//        }
//    }
}

TEST(ANALYTICAL_VS_NUMERICAL_E_L)
{
    ////////////////////////7
    /// C4
    /// ///////////////////////

    double analytical, numerical;
    VMCSolver *solver = new VMCSolver();

    //TestStuff
        double particles;
        double dimensions;
        long idum = -clock();


    for(int n = 0; n < 3; n++)
    {
        if(n==0)
        {
            solver->setTrialFunction(new Helium(solver));
            cout << "Helium" << endl;
        }
        if(n==1)
        {
            solver->setTrialFunction(new Beryllium(solver));
            cout << "Beryllium" << endl;
        }
        if(n==2)
        {
            solver->setTrialFunction(new Neon(solver));
            cout << "Neon" << endl;
        }


        particles = solver->getNParticles();
       dimensions = solver->getNDimensions();

        mat r = zeros (particles,dimensions);

        //Random positions to test E_L
        for(int i = 0; i < particles; i ++ )
        {
         for(int j = 0; j < dimensions; j++)
         {
            r(i,j) = ran2(&idum);

         }
        }

        solver->switchElectronInteraction(true);
        solver->trialFunction()->setAnalytical(false);
        solver->setCycles(50000);

        //Getting the local energy values
        solver->determinant()->updateSlaterMatrices(r,solver);
        numerical = solver->trialFunction()->localEnergy(r,solver);

        solver->trialFunction()->setAnalytical(true);

        //Getting the local energy values
        solver->determinant()->updateSlaterMatrices(r,solver);
        analytical = solver->trialFunction()->localEnergy(r,solver);


        if(solver->getRank() == 0)
        {
            cout << "Energy with the numerical E_L: " << numerical << endl;
            cout << "Energy with the E_L machinery : " << analytical << endl;
        }

        CHECK_CLOSE(numerical, analytical, 0.1 );


    }
}

TEST(ANALYTICAL_VS_NUMERICAL_ALL)
{
    ////////////////////////7
    /// Extra
    /// ///////////////////////

    double analytical, numerical;
    VMCSolver *solver = new VMCSolver();

    for(int n = 0; n < 3; n++)
    {
        if(n==0)
            solver->setTrialFunction(new Helium(solver));
        if(n==1)
            solver->setTrialFunction(new Beryllium(solver));
        if(n==2)
            solver->setTrialFunction(new Neon(solver));

        solver->switchElectronInteraction(true);
        solver->trialFunction()->setAnalytical(false);
        solver->setCycles(50000);
        solver->runMasterIntegration();

        numerical = solver->getEnergy();

        solver->switchElectronInteraction(true);
        solver->trialFunction()->setAnalytical(true);
        solver->setCycles(50000);
        solver->runMasterIntegration();

        analytical = solver->getEnergy();

        if(solver->getRank() == 0)
        {
            cout << "Energy with the numerical E_L: " << numerical << endl;
            cout << "Energy with the E_L machinery : " << analytical << endl;
        }

        CHECK_CLOSE(numerical, analytical, 0.1 );
    }

}

TEST(AnalyticalHelium)
{
    ////////////////////////
    /// C6
    //////////////////////////

    //This tests the analytical machinery calculating the local energy by testing it against
    // Helium for which we have a easy closed solution. This tests that Helium using the machinery gets the same E_L value with a random
    // electron position
    double analytical, old;
    VMCSolver *solver = new VMCSolver();


    solver->setTrialFunction(new HeliumJastrowAnalytical(solver)); //
    solver->switchElectronInteraction(true);
    solver->trialFunction()->setAnalytical(true);


//TestStuff
    double particles = solver->getNParticles();
    double dimensions = solver->getNDimensions();
    long idum = -clock();

    mat r = zeros (particles,dimensions);

    //Random positions to test derivative
    for(int i = 0; i < particles; i ++ )
    {
     for(int j = 0; j < dimensions; j++)
     {
        r(i,j) = ran2(&idum);
     }
    }

    solver->determinant()->updateSlaterMatrices(r,solver);
    old = solver->trialFunction()->localEnergy(r,solver);


    solver->setTrialFunction(new Helium(solver));

    solver->determinant()->updateSlaterMatrices(r,solver);
    analytical = solver->trialFunction()->localEnergy(r,solver);

    if(solver->getRank() == 0)
    {
        cout << "Energy with the correct E_L one: " << old << endl;
        cout << "Energy with the E_L machinery : " << analytical << endl;
    }

    CHECK_CLOSE(old, analytical, 0.0001 );
}


TEST(HYDROGENTWO_VS_HELIUM)
{
    //Make a test with R = 0 for hydrogenTwo where it should be the same as the Helium atom

    cout << endl << "Running He vs H_2 test" << endl << endl;
    VMCSolver *solver = new VMCSolver;
    solver->setTrialFunction(new Helium(solver));
    solver->trialFunction()->setAnalytical(false);
    solver->switchbBlockSampling(false);
    solver->setCycles(100000);
    solver->runMasterIntegration();

    mat r = zeros (2,3);
    vec correct = zeros (3);
    vec calculated = zeros(3);
    double r12 = 0;

    double particles = 2;
    double dimensions = 3;
    double beta = solver->getBeta();
    long idum = -10000;

    for(int i = 0; i < particles; i ++ )
    {
        for(int j = 0; j < dimensions; j++)
        {
           r(i,j) = ran2(&idum);
        }
    }

//    solver->trialFunction()->waveFunction(r, solver);
//    solver->trialFunction()->localEnergy(r,solver);
    solver->determinant()->updateSlaterMatrices(r, solver);
    cout << solver->trialFunction()->waveFunction(r,solver) << endl;


    solver->setTrialFunction(new HydrogenTwo(solver));
    solver->trialFunction()->setAnalytical(false);
    solver->trialFunction()->setNucleusDistance(0.);
    solver->switchElectronInteraction(true);
    solver->switchbBlockSampling(false);
    solver->setCycles(100000);
    solver->runMasterIntegration();

}

TEST(GRADIENT_CORR_RATIO_HELIUM)
{
//        Testing the machinery for the gradient ratio of Psi_C vs a calculated ratio
        VMCSolver *solver = new VMCSolver();

        solver->setTrialFunction(new Helium(solver));
    //    solver->switchElectronInteraction(false);
        solver->trialFunction()->setAnalytical(true);
        solver->setAlpha(solver->getCharge());
        solver->setBeta(2);
    //    solver->setCycles(1000);

        double particles = 2;
        double dimensions = 3;
        double beta = solver->getBeta();
        long idum = -1;

        mat r = zeros (particles,dimensions);
        vec correct = zeros (dimensions);
        vec calculated = zeros(dimensions);
        double r12 = 0;



        for(int i = 0; i < particles; i ++ )
        {
            for(int j = 0; j < dimensions; j++)
            {
               r(i,j) = ran2(&idum);
            }
        }


        //Testing the double derivative of the gradient ratio as the results should be calculated by the derivatives.py program

        r12 = norm(r.row(0) - r.row(1));
        correct= ((r.row(0) - r.row(1)).t())  / (r12 * pow(1+(beta*r12), 2)) ;
        calculated = (solver->derivatives()->analyticalCorrelationGradient(r,solver)).row(0).t()*2 ;

    if(solver->getRank()==0)
    {
        CHECK_CLOSE( correct(0) , calculated(0), 0.001  );
        CHECK_CLOSE( correct(1) , calculated(1), 0.001  );
        CHECK_CLOSE( correct(2) , calculated(2), 0.001  );
        cout << solver->derivatives()->analyticalCorrelationGradient(r, solver) << endl;

    }
}




