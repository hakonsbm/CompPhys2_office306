# CompPhys2_office306


The Report is in the Report folder.

Variational Monte Carlo program and studies for Helium and Beryllium

# Running the program

To run it in QT creator make sure that the run settings is set according to this guide

http://dragly.org/2012/03/14/developing-mpi-applications-in-qt-creator/


The program is run by opening the vmc.pro file with for example qt-creator.

Set run arguments to
```
-n [p] vmc [Atom] [Test] [Cycles]
```
`p` is number of processes.

`Atom` is the trialfunction. Options:
```
HeliumSimpleAnalytical   
HeliumSimpleNumerical    
HeliumJastrowAnalytical
HeliumJastrowNumerical  
Beryllium
Neon
```
`Test` is the test to be run. Options:
```
runWithDiffConstants
runSIWithDiffTimesteps
runBlockingSampledRun
runCompareAnalytical
runDiffNCycles
runFindAlphaBeta
runCompareParallelize
```
`Cycles` is the number of cycles in the Monte Carlo solver.

There is also a program called energyLevels.py that makes different plots og the various data produced by the main program.



























