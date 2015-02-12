#include "vmcsolver.h"
#include "lib.h"
#include "trialFunctions/trialfunction.h"
#include "trialFunctions/heliumsimplenumerically.h"

#include <armadillo>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>


using namespace arma;
using namespace std;

//test push

ofstream ofile;

VMCSolver::VMCSolver() :
    nDimensions(3),
    charge(2),
    stepLength(0.1),
    nParticles(2),
    h(0.001),
    h2(1000000),
    idum(-1),
    alpha(0.9*charge),
    beta(1.01),
    nCycles(1000000)
{


}

void VMCSolver::runMonteCarloIntegration()
{
    char const * outfilename = "out4-4.d";
    int acc_moves;
    double waveFunctionOld, waveFunctionNew;

    double energySum, energySquaredSum, deltaE, step_min;

    int moves, acc_moves_old;

    double beta_old = beta;


    double d_alpha = 0.1;
    double alpha_max = alpha; // 1.2*charge;
    //vec alpha_e((alpha_max-alpha)/d_alpha);
    int alpha_count = 0;
    double d_beta = 0.01;
    double beta_max = 1.5;
    //vec alpha_e((alpha_max-alpha)/d_alpha);
    //int alpha_count = 0;

    ofile.open(outfilename);
    for (alpha; alpha <= alpha_max; alpha+= d_alpha)
    {
        cout << alpha << " " << alpha_max << " " << d_alpha << endl;
        for (beta = beta_old; beta <= beta_max; beta+= d_beta)
        {
            cout << "\nAlpha/beta: " << alpha << "/" << beta << endl;

            rOld = zeros<mat>(nParticles, nDimensions);
            rNew = zeros<mat>(nParticles, nDimensions);


            waveFunctionOld = 0;
            waveFunctionNew = 0;
            energySum = 0;
            energySquaredSum = 0;
            step_min = nCycles;
            moves = 0;
            acc_moves_old = nCycles;
            stepLength = 0.1;

            // initial trial positions
            for(int i = 0; i < nParticles; i++) {
                for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
                }
            }
            rNew = rOld;

            // find optimal steplength
            for (stepLength; stepLength <= 5; stepLength += 0.1){
                acc_moves = 0;
                waveFunctionOld = 0;
                waveFunctionNew = 0;

                energySum = 0;
                energySquaredSum = 0;
                // initial trial positions
                for(int i = 0; i < nParticles; i++) {
                    for(int j = 0; j < nDimensions; j++) {
                        rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
                    }
                }
                rNew = rOld;

                for(int cycle = 0; cycle < nCycles/100; cycle++) {
                    waveFunctionOld = waveFunction(rOld);
                    for(int i = 0; i < nParticles; i++) {
                        for(int j = 0; j < nDimensions; j++) {
                            rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
                        }

                        waveFunctionNew = waveFunction(rNew);
                        if(ran2(&idum) <= (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
                            acc_moves += 1;
                            for(int j = 0; j < nDimensions; j++) {
                                rOld(i,j) = rNew(i,j);
                                waveFunctionOld = waveFunctionNew;
                            }
                        } else {
                            acc_moves -= 1;
                            for(int j = 0; j < nDimensions; j++) {
                                rNew(i,j) = rOld(i,j);
                            }
                        }
                    }
                }
                //cout << "Steplength: " << stepLength  << "  " << acc_moves << endl;
                if(abs(acc_moves) < acc_moves_old) {
                    step_min = stepLength;
                    acc_moves_old = acc_moves;
                }
            }

            stepLength = step_min;

            cout << "Steplength: " << stepLength  << "  " << acc_moves << endl;
            waveFunctionOld = 0;
            waveFunctionNew = 0;
            acc_moves = 0;
            energySum = 0;
            energySquaredSum = 0;
            // initial trial positions
            for(int i = 0; i < nParticles; i++) {
                for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
                }
            }
            rNew = rOld;
            // loop over Monte Carlo cycles
            for(int cycle = 0; cycle < nCycles; cycle++) {

                // Store the current value of the wave function
                waveFunctionOld = waveFunction(rOld);

                // New position to test
                for(int i = 0; i < nParticles; i++) {
                    for(int j = 0; j < nDimensions; j++) {
                        rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
                    }

                    // Recalculate the value of the wave function
                    waveFunctionNew = waveFunction(rNew);

                    // Check for step acceptance (if yes, update position, if no, reset position)
                    if(ran2(&idum) <= (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
                        acc_moves += 1;
                        for(int j = 0; j < nDimensions; j++) {
                            rOld(i,j) = rNew(i,j);
                            waveFunctionOld = waveFunctionNew;
                        }
                    } else {
                        for(int j = 0; j < nDimensions; j++) {
                            rNew(i,j) = rOld(i,j);
                        }
                    }
                    moves += 1;
                    // update energies
                    deltaE = localEnergy(rNew);
                    energySum += deltaE;
                    energySquaredSum += deltaE*deltaE;
                }
            }
            double energy = energySum/(nCycles * nParticles);
            double energySquared = energySquaredSum/(nCycles * nParticles);
            cout << "Energy: " << energy << " Energy (squared sum): " << energySquared << endl;
            cout << "Moves: " << moves << endl;
            cout << "Accepted moves: " << acc_moves << endl;
            cout << "Ratio: " << acc_moves/(double)moves << endl;
            ofile << setw(15) << setprecision(8) << alpha;
            ofile << setw(15) << setprecision(8) << beta;
            ofile << setw(15) << setprecision(8) << energy << endl;
            //alpha_e[alpha_count] = energy;
            alpha_count += 1;
        }
    }
    ofile.close();
}

double VMCSolver::localEnergy(const mat &r)
{
    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);

    rPlus = rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    double waveFunctionCurrent = waveFunction(r);

    // Kinetic energy

    double kineticEnergy = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = waveFunction(rMinus);
            waveFunctionPlus = waveFunction(rPlus);
            kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }
    kineticEnergy = 0.5 * h2 * kineticEnergy / waveFunctionCurrent;

    // Potential energy
    double potentialEnergy = 0;
    double rSingleParticle = 0;
    for(int i = 0; i < nParticles; i++) {
        rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j)*r(i,j);
        }
        potentialEnergy -= charge / sqrt(rSingleParticle);
    }
    // Contribution from electron-electron potential
    double r12 = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = i + 1; j < nParticles; j++) {
            r12 = 0;
            for(int k = 0; k < nDimensions; k++) {
                r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
            }
            potentialEnergy += 1 / sqrt(r12);
        }
    }

    return kineticEnergy + potentialEnergy;
}

// closed form
//double VMCSolver::localEnergy(const mat &r)
//{
//    mat rPlus = zeros<mat>(nParticles, nDimensions);
//    mat rMinus = zeros<mat>(nParticles, nDimensions);

//    rPlus = rMinus = r;

//    double waveFunctionMinus = 0;
//    double waveFunctionPlus = 0;

//    double waveFunctionCurrent = waveFunction(r);
//    double EL1, EL2;
//    double r12 = 0;
//    for(int i = 0; i < nParticles; i++) {
//        for(int j = i + 1; j < nParticles; j++) {
//            r12 = 0;
//            for(int k = 0; k < nDimensions; k++) {
//                r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
//            }
//        }
//    }
//    EL1 = (alpha - charge) * ()

//    return kineticEnergy + potentialEnergy;
//}

// Trial wavefunction T1
//double VMCSolver::waveFunction(const mat &r)
//{
//    vec rpos(nParticles);
//    for(int i = 0; i < nParticles; i++) {
//        double rSingleParticle = 0;
//        for(int j = 0; j < nDimensions; j++) {
//            rSingleParticle += r(i,j) * r(i,j);
//        }
//        rpos[i] = sqrt(rSingleParticle);
//    }
//    return exp(-accu(rpos) * alpha);
//}

// Trial wavefunction T2
 double VMCSolver::waveFunction(const mat &r)
 {
     //double r12;
     vec rpos(nParticles);
     for(int i = 0; i < nParticles; i++) {
         double rSingleParticle = 0;
         for(int j = 0; j < nDimensions; j++) {
             rSingleParticle += r(i,j) * r(i,j);
         }
         rpos[i] = sqrt(rSingleParticle);
     }
     // assuming 2 particles
     //   (ta fra elektron-elektron pot.)
     //r12 = abs(rpos[0] - rpos[1]);
     double r12 = 0;
     for(int i = 0; i < nParticles; i++) {
         for(int j = i + 1; j < nParticles; j++) {
             r12 = 0;
             for(int k = 0; k < nDimensions; k++) {
                 r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
             }
         }
     }
    

     return exp(-accu(rpos) * alpha) * exp(r12 / (2.0*(1 + beta * r12))) ;
 }
