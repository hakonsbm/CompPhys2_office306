#include "vmcsolver.h"
#include "trialFunctions/trialfunction.h"
#include "trialFunctions/heliumjastrowanalytical.h"
#include "trialFunctions/heliumjastrownumerical.h"
#include "trialFunctions/heliumsimpleanalytical.h"
#include "trialFunctions/heliumsimplenumerically.h"

#include <iostream>

using namespace std;

int main() {
    VMCSolver *solver = new VMCSolver();
    solver->setTrialFunction(new HeliumSimpleAnalytical());
    if(solver->trialFunction()->simpleFlag)


    double alpha_max = 1.2*solver->getCharge();
    double beta_max = 1.5;
    double d_alpha = 0.1;
    double d_beta = 0.01;
    cout << alpha_max << endl;

    for(double alpha = 0.9*solver->getCharge(); alpha <= 5; alpha += d_alpha) {
        solver->setAlpha(alpha);
        if(solver->trialFunction()->simpleFlag) {
            solver->calculateOptimalSteplength();
            solver->runMonteCarloIntegration();
        }
        else {
            for(double beta = 1.01; beta <= beta_max; beta += d_beta) {
                solver->setBeta(beta);
                solver->calculateOptimalSteplength();
                solver->runMonteCarloIntegration();
            }
        }
    }
    return 0;
}
