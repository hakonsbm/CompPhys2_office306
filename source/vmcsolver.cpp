#include "vmcsolver.h"
#include "lib.h"
#include "trialFunctions/trialfunction.h"
#include "trialFunctions/heliumsimplenumerical.h"

#include <armadillo>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <mpi.h>


using namespace arma;
using namespace std;
//ofstream outfile;

extern ofstream outfile;
extern ofstream samplefile;


VMCSolver::VMCSolver():
    nDimensions(3),
    charge(2),
    stepLength(0.005),
    nParticles(2),
    h(0.001),
    h2(1000000),
    idum(clock()),
    nCycles(10000),
    D(0.5),
    my_rank(0)

{
    switchElectronInteraction(true);
    initiateDerivatives(new Derivatives);
    initiateSlaterDeterminant(new SlaterDeterminant);
}



void VMCSolver::runMasterIntegration()
{


    //MPI initializations
        int numprocs;
//        MPI_Init(&m_nargs, &m_args);

        MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
//        cout << "Hello world, I have rank " << my_rank << " out of "
//        << numprocs << endl;
        double totalEnergy = 0;
        double totalEnergySquared = 0;
        double totalEnergyVar = 0;
        double totalMoves = 0;
        double totalRatio = 0;
        double totalAverageR12 = 0;


        //nCycles = nCycles/numprocs;
        runMonteCarloIntegrationIS();

        MPI_Reduce(&m_energy, &totalEnergy, 1, MPI_DOUBLE,
                      MPI_SUM, 0 , MPI_COMM_WORLD);
        MPI_Reduce(&m_energyVar, &totalEnergyVar, 1, MPI_DOUBLE,
                    MPI_SUM, 0 , MPI_COMM_WORLD);
        MPI_Reduce(&m_energySquared, &totalEnergySquared, 1, MPI_DOUBLE,
                   MPI_SUM, 0 , MPI_COMM_WORLD);
        MPI_Reduce(&m_moves, &totalMoves, 1, MPI_DOUBLE,
                   MPI_SUM, 0 , MPI_COMM_WORLD);
        MPI_Reduce(&m_ratio, &totalRatio, 1, MPI_DOUBLE,
                   MPI_SUM, 0 , MPI_COMM_WORLD);
        MPI_Reduce(&m_averageR12, &totalAverageR12, 1, MPI_DOUBLE,
                   MPI_SUM, 0 , MPI_COMM_WORLD);

        //A lot of the data should be averaged over the results from the different threads
        totalEnergy /=numprocs;
        totalEnergyVar /= numprocs;
        totalEnergySquared /= numprocs;
        totalAverageR12 /= numprocs;
        totalRatio /= numprocs;

        if(my_rank == 0)
        {
            cout << "Energy: " << totalEnergy << " Energy (squared sum): " << totalEnergySquared << endl;
            cout << "Variance: " << totalEnergyVar << endl;
            cout << "Moves: " << totalMoves << endl;
            cout << "Ratio: " <<  totalRatio << endl;
            cout << "Alpha: " << m_alpha << " and beta: " << m_beta << endl;
            cout << "Average distance between the electrons: " << totalAverageR12 << endl;
            cout << "Steplength: " << stepLength << endl;

            //Write results to file
            outfile << setw(15) << setprecision(8) << totalEnergy;
            outfile << setw(15) << setprecision(8) << totalEnergySquared;
            outfile << setw(15) << setprecision(8) << totalEnergyVar;
            outfile << setw(15) << setprecision(8) << m_alpha;
            outfile << setw(15) << setprecision(8) << m_beta;
            outfile << setw(15) << setprecision(8) << totalAverageR12;
            outfile << setw(15) << setprecision(8) << stepLength;
            outfile << setw(15) << nCycles << endl;
        }

//        // End MPI
//        MPI_Finalize ();


}



void VMCSolver::runMonteCarloIntegrationIS() {

    //Make sure that the threads have different seeds
    idum = clock();

//    char const *outfilename = "../source/outfiles/Test";// + trialFunction()->m_outfileName;
    int acc_moves = 0;
    int moves = 0;
    double waveFunctionOld = 0;
    double waveFunctionNew = 0;
    double energySum = 0;
    double energySquaredSum = 0;
    double deltaE = 0;
    double r12 = 0;
    double averageR12 = 0;

    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);
    QForceOld = zeros<mat>(nParticles, nDimensions);
    QForceNew = zeros<mat>(nParticles, nDimensions);


    //initial trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = GaussianDeviate(&idum)*sqrt(stepLength);
        }
    }

    //Set up the Slater Matrices after the move
    determinant()->updateSlaterMatrices(rOld,this);

    rNew = rOld;
    //loop over Monte Carlo cycles
    int print_cycle = 0;
    for(int cycle = 0; cycle < nCycles; cycle++) {

//        determinant()->updateSlaterMatrices(rOld,this);


        if(my_rank == 0)
        {
            if(cycle == print_cycle)
            {
                cout << (double)cycle*100./nCycles << " %" << endl;
                print_cycle +=(double) nCycles/10;
            }
        }


        //Store the current value of the wave function
        waveFunctionOld = trialFunction()->waveFunction(rOld, this);

//        QuantumForce(rOld, QForceOld);
//        cout << "Old quantumForce "<< endl << QForceOld << endl;
        derivatives()->numericalGradient(QForceOld,rOld, this);
//         cout << "New quantumForce "<< endl << QForceOld << endl;


        QForceOld = 2.*QForceOld*h/waveFunctionOld;
        //New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + GaussianDeviate(&idum)*sqrt(stepLength)+QForceOld(i,j)*stepLength*D;
            }
            //  for the other particles we need to set the position to the old position since
            //  we move only one particle at the time
            for (int k = 0; k < nParticles; k++) {
                if ( k != i) {
                    for (int j=0; j < nDimensions; j++) {
                        rNew(k,j) = rOld(k,j);
                    }
                }
            }
            //Recalculate the value of the wave function
            determinant()->updateSlaterMatrices(rNew,this);

            waveFunctionNew = trialFunction()->waveFunction(rNew, this);
            QuantumForce(rNew, QForceNew);
            QForceNew = QForceNew*h/waveFunctionNew; // possible typo





            //  we compute the log of the ratio of the greens functions to be used in the
            //  Metropolis-Hastings algorithm
            GreensFunction = 0.0;
            for (int j=0; j < nDimensions; j++) {
                GreensFunction += 0.5*(QForceOld(i,j)+QForceNew(i,j))*
                        (D*stepLength*0.5*(QForceOld(i,j)-QForceNew(i,j))-rNew(i,j)+rOld(i,j));
            }
            GreensFunction = exp(GreensFunction);

            // The Metropolis test is performed by moving one particle at the time
            MHRatio = GreensFunction*(waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld);
            if(ran2(&idum) <= MHRatio) {
                acc_moves += 1;
                for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = rNew(i,j);
                    QForceOld(i,j) = QForceNew(i,j);
                    waveFunctionOld = waveFunctionNew;
//                    determinant()->updateSlaterMatrices(rNew,this); //Updating the matrices after moving the particle :)

                }

            }
            else {
                for(int j = 0; j < nDimensions; j++) {
                   rNew(i,j) = rOld(i,j);
                   QForceNew(i,j) = QForceOld(i,j);
                   determinant()->updateSlaterMatrices(rOld,this); //Updating the matrices after moving the particle :)
                }
            }


            moves += 1;
            //update energies

            deltaE = trialFunction()->localEnergy(rNew, this);
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;

        }
        //we need to find the average value of r12
        r12 = 0;

        for(int k = 0; k < nDimensions; k++) {

            r12 += (rNew(0,k) - rNew(1,k)) * (rNew(0,k) - rNew(1,k));
        }
        averageR12 += sqrt(r12);



        if (m_blockSampling){// &&  cycle % 10 == 0) {
            samplefile << setw(15) << setprecision(8) << deltaE;
            samplefile << setw(15) << setprecision(8) << energySquaredSum;
            samplefile << setw(15) << setprecision(8) << sqrt(r12);
            samplefile << setw(15) << setprecision(8) << 1;

            for(int i = 0; i < nParticles; i++)
            {
                    samplefile << setw(15) << setprecision(8) << rNew(i,0);
                    samplefile << setw(15) << setprecision(8) << rNew(i,1);
                    samplefile << setw(15) << setprecision(8) << rNew(i,2) ;
            }
            samplefile << endl;
        }
    }

    if(m_blockSampling)
    {
        samplefile << "#Alpha: " << m_alpha << " and beta: " << m_beta << endl;
        cout << "blockSampling" << endl;
    }

    m_energy = energySum/(nCycles * nParticles);
    m_energySquared = energySquaredSum/(nCycles * nParticles);
    m_energyVar = sqrt((m_energySquared - m_energy*m_energy) / nCycles);
    m_averageR12 = averageR12 / (double) nCycles;
    m_ratio = (double) acc_moves/(double) moves;
    m_moves = moves;

}

void VMCSolver::QuantumForce(const mat &r, mat &QForce)
{
    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);


    rPlus = rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = trialFunction()->waveFunction(rMinus, this);
            waveFunctionPlus = trialFunction()->waveFunction(rPlus, this);
            QForce(i,j) =  (waveFunctionPlus-waveFunctionMinus);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }
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

void VMCSolver::runMonteCarloIntegration() {
    int acc_moves = 0;
    int moves = 0;
    double waveFunctionOld = 0;
    double waveFunctionNew = 0;
    double energySum = 0;
    double energySquaredSum = 0;
    double deltaE = 0;
    double r12 = 0;
    double averageR12 = 0;

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
    int print_cycle = 0;

    for(int cycle = 0; cycle < nCycles; cycle++) {
        if(cycle == print_cycle)
        {
            cout << (double)cycle*100./nCycles << " %" << endl;
            print_cycle += (double) nCycles/4;
        }
        //Store the current value of the wave function
        waveFunctionOld = trialFunction()->waveFunction(rOld, this);

        //New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
            }
//            determinant()->updateSlaterMatrices(rNew,this); //Updating |D| and |D|⁻1, should only be done once
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
//            determinant()->updateSlaterMatrices(rNew,this);

            deltaE = trialFunction()->localEnergy(rNew, this);
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
        }
        //we need to find the average value of r12
        r12 = 0;
//        for(int k = 0; k < nDimensions; k++) {
//            r12 += (rNew(0,k) - rNew(1,k)) * (rNew(0,k) - rNew(1,k));
//        }
        averageR12 += sqrt(r12);

        if (m_blockSampling &&  cycle % 10 == 0) {
            samplefile << setw(15) << setprecision(8) << deltaE;
            samplefile << setw(15) << setprecision(8) << deltaE*deltaE;
            samplefile << setw(15) << setprecision(8) << sqrt(r12) << endl;

            for(int i = 0; i < nParticles; i++ )
            {
                    samplefile << setw(15) << setprecision(8) << rNew(i,0);
                    samplefile << setw(15) << setprecision(8) << rNew(i,1);
                    samplefile << setw(15) << setprecision(8) << rNew(i,2);
            }
        }
    }
    if(m_blockSampling)
    {
        samplefile << "#Alpha: " << m_alpha << " and beta: " << m_beta << endl;
        cout << "blockSampling" << endl;
    }
    double energy = energySum/(nCycles * nParticles);
    double energySquared = energySquaredSum/(nCycles * nParticles);
    m_energyVar = sqrt((energySquared - energy*energy) / nCycles);  //Should we add this /(nCycles * nParticles)
    averageR12 /= (double) nCycles;

    //Storing the energy and variance calculated, used in searching for the best fit for alpha and beta
    storeEnergy(energy);
//    storeVariance(m_energyVar); //Not necessary it was already calculated

    cout << "Energy: " << energy << " Energy (squared sum): " << energySquared << endl;
    cout << "Variance: " << m_energyVar << endl;
    cout << "Moves: " << moves << endl;
    cout << "Accepted moves: " << acc_moves << endl;
    cout << "Ratio: " << (double) acc_moves/(double) moves << endl;
    cout << "Alpha: " << m_alpha << " and beta: " << m_beta << endl;
    cout << "Average distance between the electrons: " << averageR12 << endl;
    cout << "Steplength: " << stepLength << endl;



    outfile << setw(15) << setprecision(8) << energy;
    outfile << setw(15) << setprecision(8) << energySquared;
    outfile << setw(15) << setprecision(8) << m_energyVar;
    outfile << setw(15) << setprecision(8) << m_alpha;
    outfile << setw(15) << setprecision(8) << m_beta;
    outfile << setw(15) << setprecision(8) << averageR12;
    outfile << setw(15) << setprecision(8) << stepLength;
    outfile << setw(15) << nCycles << endl;
}
