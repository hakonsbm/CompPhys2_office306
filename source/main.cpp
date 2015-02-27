#include "vmcsolver.h"
#include "trialFunctions/trialfunction.h"
#include "trialFunctions/heliumjastrowanalytical.h"
#include "trialFunctions/heliumjastrownumerical.h"
#include "trialFunctions/heliumsimpleanalytical.h"
#include "trialFunctions/heliumsimplenumerical.h"
#include "trialFunctions/beryllium.h"

#include <iostream>
#include <time.h>


using namespace std;
ofstream outfile;
ofstream samplefile;

//This runs through with several different alpha and beta values and prints the results
void runWithDiffConstants(VMCSolver *solver);


int main() {
    // Choices for the trialFunctions to use are:
    //
    //HeliumSimpleAnalytical:   alpha = 1.62    beta = 0
    //HeliumSimpleNumerical:    alpha = 1.7     beta = 0
    //HeliumJastrowAnalytical:  alpha = 1.8     beta = 1.05
    //HeliumJastrowNumerical:   alpha = 1.8     beta = 1.05
    //Beryllium:                alpha = 4.0     beta = 0.31

    VMCSolver *solver = new VMCSolver();
    solver->setTrialFunction(new HeliumJastrowAnalytical()); // Beryllium

    //Opens the file that the relevant wavefunction should be written to, this file is then written to in the
    //vmcSolver class
    char const * outfilePath = (string("../source/outfiles/") + solver->trialFunction()->m_outfileName).c_str();
    char const * samplefilePath = (string("../source/outfiles/") + solver->trialFunction()->m_outfileName + string("_samples")).c_str();
    outfile.open(outfilePath);
    samplefile.open(samplefilePath);


    //Enable this if you want to calculate for all the different alpha and beta values to find the best ones.
    //Look for the program energyLevels.py to find which values weere the best
//    runWithDiffConstants(solver);

////////////////////////////////////////////////////////////////////////////////////7
//Enable the part below if you want the alph and beta values to run with (Good values below)

//////////////////////////////////////////////////////////////////////////////////////


   solver->setAlpha(1.8);
   solver->setBeta(1.05);

//  solver->calculateOptimalSteplength();
//    solver->runMonteCarloIntegration();

   solver->runMonteCarloIntegrationIS();


    cout << "\nWriting to " << outfilePath << endl;
    outfile.close();
    samplefile.close();


    return 0;
}

void runWithDiffConstants(VMCSolver *solver)
{
    double alpha_max = 1.0* solver->getCharge();
    double beta_max = 0.49;
    double d_alpha = 0.02;
    double d_beta = 0.02;

    solver->calculateOptimalSteplength();

    clock_t start, end;     //To keep track of the time


    for(double alpha = alpha = 1.0* solver->getCharge(); alpha <= alpha_max; alpha += d_alpha) {
        solver->setAlpha(alpha);
        if(solver->trialFunction()->simpleFlag) {


            start = clock();
               // solver->calculateOptimalSteplength();
            end = clock();

            //double timeOptimalStepLength = 1.0*(end - start)/CLOCKS_PER_SEC;

            start = clock();
                solver->runMonteCarloIntegrationIS();
            end = clock();

            double timeRunMonte= 1.0*(end - start)/CLOCKS_PER_SEC;

            //cout << "Time to find Optimal Steplength: " << timeOptimalStepLength << endl;
            cout << "Time to run Monte Carlo: " << timeRunMonte << endl;
        }
        else {
            for(double beta = 0.4; beta <= beta_max; beta += d_beta) {
                solver->setBeta(beta);

                start = clock();
                    //solver->calculateOptimalSteplength();
                end = clock();

                //double timeOptimalStepLength = 1.0*(end - start)/CLOCKS_PER_SEC;

                start = clock();
                    solver->runMonteCarloIntegrationIS();
                end = clock();

                double timeRunMonte= 1.0*(end - start)/CLOCKS_PER_SEC;

                //cout << "Time to find Optimal Steplength: " << timeOptimalStepLength << endl;
                cout << "Time to run Monte Carlo: " << timeRunMonte << endl;


            }
        }
    }
}
