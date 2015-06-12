
#include "vmcsolver.h"
#include "trialFunctions/trialfunction.h"
#include "trialFunctions/heliumjastrowanalytical.h"
#include "trialFunctions/heliumjastrownumerical.h"
#include "trialFunctions/heliumsimpleanalytical.h"
#include "trialFunctions/heliumsimplenumerical.h"
#include "trialFunctions/helium.h"
#include "trialFunctions/hydrogen.h"
#include "trialFunctions/beryllium.h"
#include "trialFunctions/neon.h"
#include "trialFunctions/hydrogentwo.h"
#include "slaterdeterminant.h"

#include <iostream>
#include <time.h>
#include <unittest++/UnitTest++.h>
//#include <omp.h>
#include <mpi.h>

using namespace std;
ofstream outfile;
ofstream samplefile;


//This runs through with several different alpha and beta values and prints the results
void runWithDiffConstants(VMCSolver *solver);
void runSIWithDiffTimesteps(VMCSolver *solver);
void runBlockingSampledRun(VMCSolver *solver);
void runCompareAnalytical(VMCSolver *solver);
void runDiffNCycles(VMCSolver *solver);
void runFindAlphaBeta(VMCSolver *solver);
void runCompareParallelize(VMCSolver * solver);
void runDiffBetaAndR(VMCSolver *solver);
void runNewtonsMethod(VMCSolver *solver);
int runTests(VMCSolver *solver);

int main(int nargs, char* args[])
{



    // Choices for the alpha and beta values that is set in the creation of the trialFunctions are:
    //
    //HeliumSimpleAnalytical:   alpha = 1.62    beta = 0
    //HeliumSimpleNumerical:    alpha = 1.7     beta = 0
    //HeliumJastrowAnalytical:  alpha = 1.8     beta = 1.05     1.843 0.34
    //HeliumJastrowNumerical:   alpha = 1.8     beta = 1.05
    //Beryllium:                alpha = 4.0     beta = 0.31
    //Neon:                     alpha = 10.22   beta = 0.091
    //HydrogenTwo:              alpha = 1.289   beta = 0.401
    //BerylliumTwo:             alpha = 3.725   beta = 0.246

    VMCSolver *solver = new VMCSolver();


    MPI_Init(&nargs, &args);
    solver->mpiArguments(nargs, args);

    //return  UnitTest::RunAllTests();

    // Command line argument parsing
    int my_rank, numprocs, nCycles;
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    nCycles = atoi(args[3])/numprocs;



    if((string)args[1]=="HeliumSimpleAnalytical") solver->setTrialFunction(new HeliumSimpleAnalytical(solver));
    else if((string)args[1]=="HeliumSimpleNumerical") solver->setTrialFunction(new HeliumSimpleNumerical(solver));
    else if((string)args[1]=="HeliumJastrowAnalytical") solver->setTrialFunction(new HeliumJastrowAnalytical(solver));
    else if((string)args[1]=="HeliumJastrowNumerical") solver->setTrialFunction(new HeliumJastrowNumerical(solver));
    else if((string)args[1]=="Beryllium") solver->setTrialFunction(new Beryllium(solver));
    else if((string)args[1]=="Neon") solver->setTrialFunction(new Neon(solver));
    else if((string)args[1]=="Helium") solver->setTrialFunction(new Helium(solver));
    else if((string)args[1]=="HydrogenTwo") solver->setTrialFunction(new HydrogenTwo(solver));
    else {if(my_rank==0) cout << args[1] << " is not a valid atom" << endl; exit(1);}

    if(my_rank==0)
    {
        cout<<"Atom: "<< args[1]<<endl;
        cout<<"Run: "<< args[2]<<endl;
        cout<<"Cycles: "<< args[3]<<endl;
    }

    solver->setCycles(nCycles);

    if((string)args[2]=="runWithDiffConstants") runWithDiffConstants(solver);
    else if((string)args[2]=="runSIWithDiffTimesteps") runSIWithDiffTimesteps(solver);
    else if((string)args[2]=="runBlockingSampledRun") runBlockingSampledRun(solver);
    else if((string)args[2]=="runCompareAnalytical") runCompareAnalytical(solver);
    else if((string)args[2]=="runDiffNCycles") runDiffNCycles(solver);
    else if((string)args[2]=="runFindAlphaBeta") runFindAlphaBeta(solver);
    else if((string)args[2]=="runCompareParallelize") runCompareParallelize(solver);
    else if((string)args[2]=="runTests") runTests(solver);
    else if((string)args[2]=="runDiffBetaAndR") runDiffBetaAndR(solver);
    else if((string)args[2]=="runNewtonsMethod") runNewtonsMethod(solver);
    else {if(my_rank==0) cout << args[2]  << " is not a valid runtype" << endl; exit(1);}

    // End MPI
    MPI_Finalize ();
    return 0;
}

void runDiffBetaAndR(VMCSolver *solver)
{
    //This is for use with the GTO orbitals where alpha is already given, and for two atom cores molecules where we will use the variation of the distance
    // between the nucleuses as a variational parameter along with beta.

    //Can only be run with HydrogenTwo and BerylliumTwo

    solver->trialFunction()->setNucleusDistance(1.4);

    solver->runMasterIntegration();


}

void runFindAlphaBeta(VMCSolver *solver)
{
    //Now the search for ALpha Beta should be improved by doing a more smart search, so it will first do a coarse mesh over the range
    //find the best fit for it, then make a new fine rmesh around the best fit and do a new search and so on.
    //It will aslo make a better fit by taking into account both the variance (which should be zero) and the energy (which should be as low as possible),
    //when the lowest energy and variance does not agree anymore the limit for the resulition of the search has been reached, if not specified in some other way.

    double alphaMin = 0.7* solver->getCharge();
    double alphaMax = 1.2* solver->getCharge();
    double betaMin = 0.3;
    double betaMax = 0.4;

    int nSteps = 3;    //Coarseness of mesh
    int nMeshes = 4;    //Number of times it should decrease the mesh
    double meshRangeAlpha = alphaMax - alphaMin;   //Used to recalculate the mesh
    double meshRangeBeta = betaMax - betaMin;

    double bestAlphaEnergy = 1000;    //The currently best values for alpha and beta, both found by minimizing energy and variance respectively
    double bestBetaEnergy = 1000;
    double bestAlphaVariance = 1000;
    double bestBetaVariance = 1000;
    double lowestVariance = 1000;
    double lowestEnergy = 1000;

    double dAlpha ;
    double dBeta ;

    bool ImportanceSampling = true;    //Set to true if you want to run with importance sampling
    solver->switchbBlockSampling(false);
//    solver->setCycles(1000000);

    //Temporary for the HydrogenTwo function
    solver->trialFunction()->setNucleusDistance(1.4);

    //Opens the file that the relevant wavefunction should be written to, this file is then written to in the
    //vmcSolver class
    string pathString = "../source/outfiles/" +  solver->trialFunction()->m_outfileName;
    char const * outfilePath = (pathString + string("_alpha_beta")).c_str();
    outfile.open(outfilePath);

    clock_t start, end;     //To keep track of the time
    for(int i = 0; i < nMeshes; i++)
    {
        //Recalculates the mesh to be finer
        dAlpha = (alphaMax-alphaMin)/ (double) nSteps;
        dBeta = (betaMax-betaMin)/ (double) nSteps;

        for(double alpha = alphaMin ; alpha <= alphaMax; alpha += dAlpha) {
            solver->setAlpha(alpha);
            if(solver->trialFunction()->simpleFlag) {
                if(ImportanceSampling)
                {
                    start = clock();
                        solver->runMasterIntegration();
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

                if( solver->getEnergyVar() < lowestVariance )
                {
                    lowestVariance = solver->getEnergyVar();
                    bestAlphaVariance = alpha;
                }

                if(solver->getEnergy() < lowestEnergy)
                {
                    lowestEnergy = solver->getEnergy();
                    bestAlphaEnergy = alpha;
                }

            } else
            {
                for(double beta = betaMin ; beta <= betaMax; beta += dBeta) {


                    solver->setBeta(beta);
                    if(ImportanceSampling)
                    {
                        start = clock();
                            solver->runMasterIntegration();
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

                    if( solver->getEnergyVar() < lowestVariance )
                    {
                        lowestVariance = solver->getEnergyVar();
                        bestAlphaVariance = alpha;
                        bestBetaVariance = beta;
                    }

                    if(solver->getEnergy() < lowestEnergy)
                    {
                        lowestEnergy = solver->getEnergy();
                        bestAlphaEnergy = alpha;
                        bestBetaEnergy = beta;
                    }
                }
            }
        }
        cout << " The lowest variance, " << lowestVariance << ", is found with alpha  " << bestAlphaVariance  << " and beta " << bestBetaVariance << endl;
        cout << " The lowest energy, " << lowestEnergy << ", is found with alpha  " << bestAlphaEnergy  << " and beta " << bestBetaEnergy << endl;

        if(bestAlphaVariance != bestAlphaEnergy)
            break;

//        if(bestBetaEnergy != bestBetaVariance)        //Not canceling the loop because of beta, since beta has a very small influence on the values and could fail because of randomness
//            break;


        //Sets upa new more finegrained mesh
        meshRangeAlpha = meshRangeAlpha/2.;
        alphaMin = bestAlphaEnergy - meshRangeAlpha/2.;
        alphaMax = bestAlphaEnergy + meshRangeAlpha/2.;


        meshRangeBeta = meshRangeBeta/2.;
        betaMin = bestBetaEnergy - meshRangeBeta/2.;
        betaMax = bestBetaEnergy + meshRangeBeta/2.;

        cout << "New mesh is; " << alphaMax << " - " << alphaMin << endl;
        cout << "and for beta:" << betaMax << " - " << betaMin << endl;


    }
}

void runWithDiffConstants(VMCSolver *solver)
{
    //Settings for which values it should be cycled over and if we want to use importance sampling or now
    int my_rank;
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    double alpha_min = 9.5;//0.7*solver->getCharge();
    double alpha_max = 10.5;//1.5* solver->getCharge();

    int nSteps = 10;

    double beta_min = 0.05;
    double beta_max = 0.13;
    double d_alpha = (alpha_max-alpha_min)/ (double) nSteps;
    double d_beta = (beta_max-beta_min)/ (double) nSteps;

    solver->switchElectronInteraction(true);
    solver->trialFunction()->setAnalytical(false);

    bool ImportanceSampling = true;    //Set to true if you want to run with importance sampling
    solver->switchbBlockSampling(false);
//    solver->setCycles(1000000);

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
                    solver->runMasterIntegration();
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
                    solver->runMasterIntegration();
                    end = clock();

                    double timeRunMonte= 1.0*(end - start)/CLOCKS_PER_SEC;
                    if(my_rank == 0) cout << "Time to run Monte Carlo: " << timeRunMonte << endl;
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
//    solver->setCycles(1000000);

    int nSteps = 100;
    double time_min = 0.01;
    double time_max = 1.;
    double dt = (time_max-time_min)/ (double) nSteps;

    double timeStep;

    solver->setAlpha(4);
    solver->setBeta(0.31);

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
    solver->switchElectronInteraction(true);
    solver->trialFunction()->setAnalytical(false);

    string pathString = "../source/outfiles/" +  solver->trialFunction()->m_outfileName;

    char const * outfilePath = (pathString + string("_blockingSamples")).c_str();

//    solver->setCycles(1000000);

    samplefile.open(outfilePath);

    solver->runMasterIntegration();

    samplefile.close();

}

void runCompareAnalytical(VMCSolver *solver)
{
    double timeRunAnalytic;
    double timeRunNumerical;


    solver->switchbBlockSampling(false);
//    solver->setCycles(10000000);

    clock_t start, end;     //To keep track of the time

    start = clock();
        solver->runMonteCarloIntegration();
    end = clock();

    timeRunAnalytic = 1.0*(end - start)/CLOCKS_PER_SEC;

    solver->setTrialFunction(new HeliumJastrowNumerical(solver));

    solver->trialFunction();

    start = clock();
        solver->runMonteCarloIntegration();
    end = clock();

    timeRunNumerical = 1.0*(end - start)/CLOCKS_PER_SEC;

    cout << "Time to calculate analytic vs numerical " << timeRunAnalytic << " vs " << timeRunNumerical << endl;
    cout << "Time run gain "  <<  (timeRunAnalytic - timeRunNumerical) / timeRunNumerical << endl;
    cout << "Time ratia " << timeRunNumerical/timeRunAnalytic << endl;

}

void runDiffNCycles(VMCSolver *solver)
{
    string pathString = "../source/outfiles/" +  solver->trialFunction()->m_outfileName;

    char const * outfilePath = (pathString + string("_nCycles")).c_str();

    int wantedCycles = 10000000;
    int nSteps = 100;
    double deltaCycles = (double) wantedCycles/nSteps;
    int nCycles = 100;

    outfile.open(outfilePath);

    for(nCycles = 100; nCycles < wantedCycles ; nCycles += deltaCycles )
    {
        solver->setCycles(nCycles);

        solver->runMonteCarloIntegrationIS();
    }

    outfile.close();

}

void runCompareParallelize(VMCSolver * solver)
{
    //TestSettings
    solver->switchElectronInteraction(true);
    solver->trialFunction()->setAnalytical(true);
//    solver->setAlpha(solver->getCharge());



    double start, end;
    int numprocs;
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    //Need to make a python script to run it with different number of nodes
    //Or we can do it manually with changing the projoect build
       start = MPI_Wtime();
       solver->runMasterIntegration();
//       solver->runMonteCarloIntegrationIS();
       end = MPI_Wtime();

    if (solver->getMy_Rank()==0)
    {
        cout << "Time used with "<< numprocs << " processes: " << end - start << endl;
    }


    return;
}

void runNewtonsMethod(VMCSolver *solver)
{

    solver->trialFunction()->setAnalytical(true);
    solver->trialFunction()->setConjugate(true);

    solver->switchbBlockSampling(false);    //This also samples the energies at each cycle to do blocking analysis on the data

    //Opens the file that the relevant wavefunction should be written to, this file is then written to in the
    //vmcSolver class
    string pathString = "../source/outfiles/" +  solver->trialFunction()->m_outfileName;

    char const * outfilePath = (pathString + string("_conjugate")).c_str();

    outfile.open(outfilePath);

    //Assuming that alpha is known, or GTO trial functions beta is what we want to minimize E_L against. So we guess a value and goes from there
    int steps = 0;
    double lowerEnd = 0;
    double lowerDerivative = 0;
    double midPoint = 0.5;
    double midDerivative = 0;
    double higherEnd = 1;
    double higherDerivative = 0;

    solver->setBeta(lowerEnd);
    solver->runMasterIntegration();
    lowerDerivative = solver->getEnergyDerivative();
    solver->setBeta(higherEnd);
    solver->runMasterIntegration();
    higherDerivative = solver->getEnergyDerivative();
    if(higherDerivative/abs(higherDerivative) == -lowerDerivative/abs(lowerDerivative)) {
        solver->setBeta(midPoint);
        while(((higherEnd-lowerEnd)*0.5 > 0.01) && (steps <= 100)) {
            solver->runMasterIntegration();
            midDerivative = solver->getEnergyDerivative();
            if(midDerivative/abs(midDerivative) == lowerDerivative/abs(lowerDerivative)) {
                lowerEnd = midPoint;
                midPoint += (higherEnd-lowerEnd)*0.5;
            }
            else if(midDerivative/abs(midDerivative) == higherDerivative/abs(higherDerivative)) {
                higherEnd = midPoint;
                midPoint -= (higherEnd-lowerEnd)*0.5;
            }
            solver->setBeta(midPoint);
            steps++;
        }
    }
    else {
        cout << "No roots can be guaranteed to be found in the interval between " << lowerEnd << " and " << higherEnd << endl;
    }




//    double beta = 1.11; //This is the initial seed
//    double h= 0.001;
//    double firstDerivative = 0;

//    solver->setBeta(beta);
//    while((abs(solver->getEnergyDerivative()) > 0.01) && (steps <= 100)) {
//        solver->runMasterIntegration();
//        firstDerivative = solver->getEnergyDerivative();
//        solver->setBeta(beta+h);
//        solver->runMasterIntegration();
//        beta -= (firstDerivative*h)/(solver->getEnergyDerivative()-firstDerivative);
//        solver->setBeta(beta);
//        cout << "D=" << firstDerivative << endl;
//        steps++;
//   }


    if(solver->getRank() == 0)
    {
        cout << "\nWriting to " << outfilePath << endl;
    }

    outfile.close();

return;
}

int runTests(VMCSolver *solver)
{
    return  UnitTest::RunAllTests();
}
