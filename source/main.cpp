#include "vmcsolver.h"
#include "trialFunctions/trialfunction.h"
#include "trialFunctions/heliumjastrowanalytical.h"
#include "trialFunctions/heliumjastrownumerical.h"
#include "trialFunctions/heliumsimpleanalytical.h"
#include "trialFunctions/heliumsimplenumerical.h"

#include <iostream>
#include <time.h>


using namespace std;
ofstream outfile;

int main() {
    VMCSolver *solver = new VMCSolver();
    solver->setTrialFunction(new HeliumSimpleAnalytical());

    //Opens the file that the relevant wavefunction should be written to, this file is then written to in the
    //vmcSolver class
    char const * outfilePath = (string("../source/outfiles/") + solver->trialFunction()->m_outfileName).c_str();
    outfile.open(outfilePath);


    double alpha_max = 1.0*solver->getCharge();
    double beta_max = 1.5;
    double d_alpha = 0.1;
    double d_beta = 0.01;

    clock_t start, end;     //To keep track of the time


    for(double alpha = 0.5*solver->getCharge(); alpha <= alpha_max; alpha += d_alpha) {
        solver->setAlpha(alpha);
        if(solver->trialFunction()->simpleFlag) {

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
        else {
            for(double beta = 1.01; beta <= beta_max; beta += d_beta) {
                solver->setBeta(beta);

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

    outfile.close();


    return 0;
}
