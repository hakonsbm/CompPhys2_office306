#include "vmcsolver.h"
#include "trialFunctions/trialfunction.h"
#include "trialFunctions/heliumjastrowanalytical.h"
#include "trialFunctions/heliumjastrownumerical.h"
#include "trialFunctions/heliumsimpleanalytical.h"
#include "trialFunctions/heliumsimplenumerical.h"
#include "trialFunctions/hydrogen.h"
#include "trialFunctions/beryllium.h"

#include <iostream>
#include <time.h>
#include <unittest++/UnitTest++.h>
#include <omp.h>


using namespace std;
ofstream outfile;
ofstream samplefile;

//This runs through with several different alpha and beta values and prints the results
void runWithDiffConstants(VMCSolver *solver);
void runSIWithDiffTimesteps(VMCSolver *solver);
void runBlockingSampledRun(VMCSolver *solver) ;




int main() {
    // Choices for the alpha and beta values that is set in the creation of the trialFunctions are:
    //
    //HeliumSimpleAnalytical:   alpha = 1.62    beta = 0
    //HeliumSimpleNumerical:    alpha = 1.7     beta = 0
    //HeliumJastrowAnalytical:  alpha = 1.8     beta = 1.05     1.843 0.34
    //HeliumJastrowNumerical:   alpha = 1.8     beta = 1.05
    //Beryllium:                alpha = 4.0     beta = 0.31

    VMCSolver *solver = new VMCSolver();
    solver->setTrialFunction(new Beryllium(solver));


    //Enable this if you want to calculate for all the different alpha and beta values to find the best ones.
    //Look for the program energyLevels.py to find which values were the best
    runWithDiffConstants(solver);

//    runSIWithDiffTimesteps(solver);
//    runBlockingSampledRun(solver) ;


//  return  UnitTest::RunAllTests();
    return 0;
}

void runWithDiffConstants(VMCSolver *solver)
{
    //Settings for which values it should be cycled over and if we want to use importance sampling or now

    double alpha_min = 0.8*solver->getCharge();
    double alpha_max = 1.2* solver->getCharge();

    int nSteps = 10;

    double beta_min = 0.3;
    double beta_max = 0.4;
    double d_alpha = (alpha_max-alpha_min)/ (double) nSteps;
    double d_beta = (beta_max-beta_min)/ (double) nSteps;

    bool ImportanceSampling = true;    //Set to true if you want to run with importance sampling
    solver->switchbBlockSampling(false);
    solver->setCycles(1000000);

    //Opens the file that the relevant wavefunction should be written to, this file is then written to in the
    //vmcSolver class
    string pathString = "../source/outfiles/" +  solver->trialFunction()->m_outfileName;

    char const * outfilePath = (pathString + string("_alpha_beta")).c_str();


    outfile.open(outfilePath);
//    samplefile.open(samplefilePath);



    clock_t start, end;     //To keep track of the time

    for(double alpha = alpha_min ; alpha <= alpha_max; alpha += d_alpha) {
        solver->setAlpha(alpha);
        if(solver->trialFunction()->simpleFlag) {
            if(ImportanceSampling)
            {
                start = clock();
                    solver->runMonteCarloIntegrationIS();
                end = clock();

                double timeRunMonte= 1.0*(end - start)/CLOCKS_PER_SEC;
                cout << "Time to run Monte Carlo: " << timeRunMonte << endl;
            }
            else {
                start = clock();
                    solver->calculateOptimalSteplength();
//                      solver->setStepLength(1.4);
                end = clock();

                double timeOptimalStepLength = 1.0*(end - start)/CLOCKS_PER_SEC;

                start = clock();
                    solver->runMonteCarloIntegration();
                end = clock();

                double timeRunMonte= 1.0*(end - start)/CLOCKS_PER_SEC;

                cout << "Time to find Optimal Steplength: " << timeOptimalStepLength << endl;
                cout << "Time to run Monte Carlo: " << timeRunMonte << endl;
            }
        }
        else {
            for(double beta = beta_min ; beta <= beta_max; beta += d_beta) {
                solver->setBeta(beta);
                if(ImportanceSampling)
                {
                    start = clock();
                        solver->runMonteCarloIntegrationIS();
                    end = clock();

                    double timeRunMonte= 1.0*(end - start)/CLOCKS_PER_SEC;
                    cout << "Time to run Monte Carlo: " << timeRunMonte << endl;
                }
                else {
                    start = clock();
                        solver->calculateOptimalSteplength();
                    end = clock();

                    double timeOptimalStepLength = 1.0*(end - start)/CLOCKS_PER_SEC;

                    start = clock();
                        solver->runMonteCarloIntegration();
                    end = clock();

                    double timeRunMonte= 1.0*(end - start)/CLOCKS_PER_SEC;

                    cout << "Time to find Optimal Steplength: " << timeOptimalStepLength << endl;
                    cout << "Time to run Monte Carlo: " << timeRunMonte << endl;
                }
            }
        }
    }

    cout << "\nWriting to " << outfilePath << endl;
    outfile.close();
//    samplefile.close();
}

void runSIWithDiffTimesteps(VMCSolver *solver)
{
    solver->switchbBlockSampling(false);
    solver->setCycles(100000);

    int nSteps = 100;
    double time_min = 0.001;
    double time_max = 1.;
    double dt = (time_max-time_min)/ (double) nSteps;

    double timeStep;

    solver->setAlpha(1.843 );
    solver->setBeta(0.34);

    solver->switchbBlockSampling(false);    //This also samples the energies at each cycle to do blocking analysis on the data

    //Opens the file that the relevant wavefunction should be written to, this file is then written to in the
    //vmcSolver class
    string pathString = "../source/outfiles/" +  solver->trialFunction()->m_outfileName;

    char const * outfilePath = (pathString + string("_timeStep")).c_str();

    outfile.open(outfilePath);

    clock_t start, end;     //To keep track of the time


    for(timeStep = time_min ; timeStep < time_max ; timeStep += dt )
    {

        solver->setStepLength(timeStep);

        start = clock();
            solver->runMonteCarloIntegrationIS();
        end = clock();

        double timeRunMonte= 1.0*(end - start)/CLOCKS_PER_SEC;

        cout << "Time to run Monte Carlo: " << timeRunMonte << endl;

    }

    outfile.close();
    cout << "\nWriting to " << outfilePath << endl;

}

void runBlockingSampledRun(VMCSolver *solver)
{
    solver->switchbBlockSampling(true);

    string pathString = "../source/outfiles/" +  solver->trialFunction()->m_outfileName;

    char const * outfilePath = (pathString + string("_blockingSamples")).c_str();

    solver->setCycles(1000000);

    samplefile.open(outfilePath);

    solver->runMonteCarloIntegrationIS();

    samplefile.close();

}

TEST(Hydrogen) {



    VMCSolver *solver = new VMCSolver();
    solver->setTrialFunction(new Hydrogen(solver));
    solver->runMonteCarloIntegration();
    CHECK_EQUAL(0., solver->getEnergyVar());

}
