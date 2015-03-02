# CompPhys2_office306

########################################
The Report is in the Report folder.

Variational Monte Carlo program and studies for Helium and Beryllium

The program is run by opening the vmc.pro file with for example qt-creator.

Then in the program there several trialfunction built in that all can be set in the main function.

    HeliumSimpleAnalytical:   
    HeliumSimpleNumerical:    
    HeliumJastrowAnalytical:  
    HeliumJastrowNumerical:  
    Beryllium:               

These are set by using the solver->setTrialFunction(new TrialFunction(solver));

Then there are several option for which test that should be run:

void runWithDiffConstants(VMCSolver *solver);
void runSIWithDiffTimesteps(VMCSolver *solver);
void runBlockingSampledRun(VMCSolver *solver);
void runCompareAnalytical(VMCSolver *solver);
void runDiffNCycles(VMCSolver *solver);

There is also a program called energyLevels.py that makes different plots og the various data produced by the main program.



























