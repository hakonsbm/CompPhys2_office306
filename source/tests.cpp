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

//    exit(0);


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

//    exit(0);

}

TEST(HeliumJastrow)
{
//    cout << "Testing The analytical version of the Helium wavefunction" << endl;

//    double oldVersion, newVersion;


//    //Testing the machinery for the gradient ratio of Psi_C vs a calcualted ratio
//    VMCSolver *solver = new VMCSolver();

//    solver->setTrialFunction(new Helium(solver));
////    solver->switchElectronInteraction(false);
////    solver->trialFunction()->setAnalytical(true);
//    solver->setAlpha(solver->getCharge());
//    solver->setBeta(2);
////    solver->setCycles(1000);

//    mat r = zeros (2,3);
//    vec correct = zeros (3);
//    vec calculated = zeros(3);
//    double r12 = 0;

//    double particles = 2;
//    double dimensions = 3;
//    double beta = solver->getBeta();
//    long idum = -1;

//    for(int i = 0; i < particles; i ++ )
//    {
//        for(int j = 0; j < dimensions; j++)
//        {
//           r(i,j) = ran2(&idum);
//        }
//    }

//    r12 = norm(r.row(0) - r.row(1));
//    correct= ((r.row(0) - r.row(1)).t())  / (r12 * pow(1+(beta*r12), 2)) ;
//    calculated = solver->derivatives()->analyticalCorrelationDerivative(r,solver) ;

//    CHECK_CLOSE( correct(0) , calculated(0), 0.001  );
//    CHECK_CLOSE( correct(1) , calculated(1), 0.001  );
//    CHECK_CLOSE( correct(2) , calculated(2), 0.001  );


//    //Checking the laplacian ratio of Psi_C vs a calculated one.
////    cout << solver->derivatives()->analyticalCorrelationDoubleDerivative(r, solver) << endl;


////    exit(0);
}

TEST(FirstOrderDerivative)
{
    VMCSolver *solver = new VMCSolver();

    solver->setTrialFunction(new Helium(solver));


    double particles = solver->getNParticles();
    double dimensions = solver->getNDimensions();
    long idum = -1;

    mat r = zeros (particles,dimensions);
    mat gradientNumerical = zeros(particles, dimensions);
    mat gradientAnalytical = zeros(particles, dimensions);


    //Random positions to test derivative
    for(int i = 0; i < particles; i ++ )
    {
     for(int j = 0; j < dimensions; j++)
     {
        r(i,j) = ran2(&idum);
     }
    }
    solver->trialFunction()->setAnalytical(false);
    solver->determinant()->updateSlaterMatrices(r,solver);
    solver->derivatives()->numericalGradient(gradientNumerical, r, solver);

    solver->trialFunction()->setAnalytical(true);
    solver->determinant()->updateSlaterMatrices(r,solver);
    solver->derivatives()->analyticalGradient(gradientAnalytical, r , solver);

    if(solver->getRank() == 0)
    {
        cout << "End results" << endl;
        cout << gradientNumerical << endl;
        cout << gradientAnalytical << endl;
    }

}

TEST(AnalyticalHelium)
{

////        Testing the machinery for the gradient ratio of Psi_C vs a calculated ratio
//        VMCSolver *solver = new VMCSolver();

//        solver->setTrialFunction(new Helium(solver));
//    //    solver->switchElectronInteraction(false);
//        solver->trialFunction()->setAnalytical(true);
//        solver->setAlpha(solver->getCharge());
//        solver->setBeta(2);
//    //    solver->setCycles(1000);

//        mat r = zeros (2,3);
//        vec correct = zeros (3);
//        vec calculated = zeros(3);
//        double r12 = 0;

//        double particles = 2;
//        double dimensions = 3;
//        double beta = solver->getBeta();
//        long idum = -1;

//        for(int i = 0; i < particles; i ++ )
//        {
//            for(int j = 0; j < dimensions; j++)
//            {
//               r(i,j) = ran2(&idum);
//            }
//        }


//        //Testing the double derivative of the gradient ratio as the results should be calculated by the derivatives.py program

//        r12 = norm(r.row(0) - r.row(1));
//        correct= ((r.row(0) - r.row(1)).t())  / (r12 * pow(1+(beta*r12), 2)) ;
//        calculated = solver->derivatives()->analyticalCorrelationDerivative(r,solver) ;
//    if(solver->getRank()==0)
//    {
//        CHECK_CLOSE( correct(0) , calculated(0), 0.001  );
//        CHECK_CLOSE( correct(1) , calculated(1), 0.001  );
//        CHECK_CLOSE( correct(2) , calculated(2), 0.001  );
//        cout << solver->derivatives()->analyticalCorrelationDerivative(r, solver) << endl;

//    }

    //Testing the laplacian ratio


////    cout << "Testing The analytical version of the Helium wavefunction" << endl;

//    double analytical, numerical, analyticalVar , numericalVar;
//    VMCSolver *solver = new VMCSolver();


////    cout << endl << "Running Helium test" << endl << endl;
//    solver->setTrialFunction(new HeliumJastrowAnalytical(solver));
//    solver->switchElectronInteraction(true);
//    solver->trialFunction()->setAnalytical(true);
//    solver->switchbBlockSampling(false);
//    solver->setCycles(2);
//    solver->runMasterIntegration();

//    numerical = solver->getEnergy();
//    numericalVar = solver->getEnergyVar();


//    solver->setTrialFunction(new Helium(solver));
//    solver->switchElectronInteraction(true);
//    solver->trialFunction()->setAnalytical(false);
//    solver->switchbBlockSampling(false);
//    solver->setCycles(2);
//    solver->runMasterIntegration();

//    analytical = solver->getEnergy();
//    analyticalVar = solver->getEnergyVar();


//    if(solver->getRank() == 0)
//    {
//        cout << "Energy with the old one: " << numerical << " and " << numericalVar << endl;
//        cout << "Energy with the new one: " << analytical << " and " << analyticalVar << endl;
//    }


}

