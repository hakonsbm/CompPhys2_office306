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

    solver->setTrialFunction(new HeliumSimpleNumerically());

    solver->calculateOptimalSteplength();
    double alpha_max = 1.2*solver->getCharge();
    double beta_max = 1.5;
    double d_alpha = 0.1;
    double d_beta = 0.01;
    for(double alpha = 0.9*solver->getCharge(); alpha <= alpha_max; alpha += d_alpha) {
        for(double beta = 1.01; beta <= beta_max; beta += d_beta) {
            solver->runMonteCarloIntegration(alpha, beta);
        }
    }
    return 0;
}
