#include "vmcsolver.h"

#include <iostream>

using namespace std;

int main()
{
    VMCSolver *solver = new VMCSolver();

    solver->setTrialFunction(new HeliumSimpleNumerically());

    solver->runMonteCarloIntegration();
    return 0;

}
