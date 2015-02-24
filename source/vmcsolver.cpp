#include "vmcsolver.h"
#include "lib.h"
#include "trialFunctions/trialfunction.h"
#include "trialFunctions/heliumsimplenumerical.h"

#include <armadillo>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>


using namespace arma;
using namespace std;
//ofstream outfile;

extern ofstream outfile;



VMCSolver::VMCSolver():
    nDimensions(3),
    charge(2),
    stepLength(0.1),
    nParticles(2),
    h(0.001),
    h2(1000000),
    idum(-1),
    nCycles(1000000)
{

}

void VMCSolver::runMonteCarloIntegration() {
//    char const *outfilename = "../source/outfiles/Test";// + trialFunction()->m_outfileName;
    int acc_moves = 0;
    int moves = 0;
    double waveFunctionOld = 0;
    double waveFunctionNew = 0;
    double energySum = 0;
    double energySquaredSum = 0;
    double deltaE = 0;

    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);

    //initial trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
        }
    }

    rNew = rOld;
    //loop over Monte Carlo cycles
    for(int cycle = 0; cycle < nCycles; cycle++) {

        //Store the current value of the wave function
        waveFunctionOld = trialFunction()->waveFunction(rOld, this);


        //New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
            }
            //Recalculate the value of the wave function
            waveFunctionNew = trialFunction()->waveFunction(rNew, this);

            //Check for step acceptance (if yes, update position, if no, reset position)
            if(ran2(&idum) <= (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
                acc_moves += 1;
                for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = rNew(i,j);
                    waveFunctionOld = waveFunctionNew;
                }
            }
            else {
                for(int j = 0; j < nDimensions; j++) {
                   rNew(i,j) = rOld(i,j);
                }
            }
            moves += 1;
            //update energies
            deltaE = trialFunction()->localEnergy(rNew, this);
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
        }
    }
    double energy = energySum/(nCycles * nParticles);
    double energySquared = energySquaredSum/(nCycles * nParticles);
    double energyVar = energySquared - energy*energy;


    cout << "Energy: " << energy << " Energy (squared sum): " << energySquared << endl;
    cout << "Variance: " << energyVar << endl;
    cout << "Moves: " << moves << endl;
    cout << "Accepted moves: " << acc_moves << endl;
    cout << "Ratio: " << (double) acc_moves/(double) moves << endl;
    cout << "Alpha: " << m_alpha << " and beta: " << m_beta << endl;

    outfile << setw(15) << setprecision(8) << energy;
    outfile << setw(15) << setprecision(8) << energySquared;
    outfile << setw(15) << setprecision(8) << m_alpha;
    outfile << setw(15) << setprecision(8) <<  m_beta << endl;
}



void VMCSolver::calculateOptimalSteplength() {
    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);

    double waveFunctionOld = 0;
    double waveFunctionNew = 0;
    double step_min = 1;
    double ratio = 0;
    double old_ratio = 1;
    int moves = 0;
    int acc_moves = 0;

    // initial trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = stepLength*(ran2(&idum) - 0.5);
        }
    }
    rNew = rOld;

    // find optimal steplength
    for (stepLength = 0.1; stepLength <= 20.0; stepLength += 0.2){
        moves = 0;
        acc_moves = 0;
        waveFunctionOld = 0;
        waveFunctionNew = 0;

        // initial trial positions
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rOld(i,j) = stepLength*(ran2(&idum) - 0.5);
            }
        }
        rNew = rOld;
        for(int cycle = 0; cycle < 1000; cycle++) {
            waveFunctionOld = trialFunction()->waveFunction(rOld, this);
            for(int i = 0; i < nParticles; i++) {
                for(int j = 0; j < nDimensions; j++) {
                    rNew(i,j) = rOld(i,j)+stepLength*(ran2(&idum) - 0.5);
                }
                waveFunctionNew = trialFunction()->waveFunction(rNew, this);
                if(ran2(&idum) <= (waveFunctionNew*waveFunctionNew)/(waveFunctionOld*waveFunctionOld)) {
                    acc_moves += 1;
                    for(int j = 0; j < nDimensions; j++) {
                        rOld(i,j) = rNew(i,j);
                        waveFunctionOld = waveFunctionNew;
                    }
                } else {
                    //acc_moves -= 1;
                    for(int j = 0; j < nDimensions; j++) {
                        rNew(i,j) = rOld(i,j);
                    }
                }
                moves += 1;
            }
        }
        ratio = (double)acc_moves/(double)moves;
        if(abs(0.5-ratio) < abs(0.5-old_ratio)) {
            step_min = stepLength;
            old_ratio = ratio;
        }
    }
    stepLength = step_min;
    cout << endl << "##############################################################################" << endl << endl;
    cout << endl << "Steplength: " << stepLength << endl;
}


