#include "berylliumtwo.h"


BerylliumTwo::BerylliumTwo(VMCSolver *solver)
{
    simpleFlag = false;
    m_analytical = false;
    m_molecule = true;

    m_outfileName = "BerylliumTwo";

//    cout << "got here" << endl;


//    solver->trialFunction()->setNucleusDistance(1.4);

//    cout << "got here" << endl;


    solver->setCharge(4);
    solver->setNParticles(8);
    solver->setAlpha(3.725);
    solver->setBeta(0.246);

    spin << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1;

}



double BerylliumTwo::waveFunction(const mat &r, VMCSolver *solver)
{
    double SD, rij, a;     //rip1 is the distance from electron i to nucleus 1
    double product = 1.0;
    double alpha = solver -> getAlpha();
    double beta = solver -> getBeta();
    double nParticles = solver->getNParticles();
    double nDimensions = solver->getNDimensions();


    //Calculate the Jastrow factor
    if(solver->getElectronInteration())
    {
        for(int i = 0; i < nParticles; i++) {
            for(int j = i + 1; j < nParticles; j++) {
                rij = 0;
                for(int k = 0; k < nDimensions; k++) {
                    rij += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
                }
                rij = sqrt(rij);
                a = spinFactor(i,j);
                product = product * exp(a*rij/(1+beta*rij));
            }
        }
    }

    SD = solver->determinant()->calculateDeterminant(r,alpha,solver); //SlaterDeterminant(r, alpha, solver);

    return SD*product;
}

double BerylliumTwo::localEnergy(const mat &r, VMCSolver *solver)
{
    double kineticEnergy, potentialEnergy;
    int nDimensions = solver->getNDimensions();
    double charge = solver->getCharge();
    int nParticles = solver->getNParticles();

    vec nucleusHalfDistance = zeros(nDimensions);
    nucleusHalfDistance(nDimensions - 1) = solver->trialFunction()->getNucleusDistance()/2.;


//    double r1p1 = norm(r.row(0) + nucleusHalfDistance.t());
//    double r1p2 = norm(r.row(0) - nucleusHalfDistance.t());
//    double r2p1 = norm(r.row(1) + nucleusHalfDistance.t());
//    double r2p2 = norm(r.row(1) - nucleusHalfDistance.t());

//    double r12 = norm(r.row(0) - r.row(1));



    if(m_analytical)
    {
    //Calculates the kinetic energy as the ratios of -1/2* ( d²/dx²|D| /|D| + 2 * (d/dx |D|/|D|)*d/dx Psi_C/Psi_C + d²/dx² Psi_C /Psi_C )
            kineticEnergy += solver->determinant()->laplacianSlaterDeterminant(r, solver);
        if(solver->getElectronInteration())
        {
            kineticEnergy += solver->derivatives()->analyticalCorrelationDoubleDerivative(r,solver);
//            cout << solver->derivatives()->analyticalCorrelationDoubleDerivative(r,solver) << endl;
//            kineticEnergy += 2*(dot(gradientSlater, gradientJastrow ));
        }

            kineticEnergy *= -1./2.;
    }
    else
        kineticEnergy =  solver->derivatives()->numericalDoubleDerivative(r, solver) / (2.*waveFunction(r, solver));


    //Time to calculate the potential energy, this is divided into three parts, potentialenergy between the nuclei,
    //potentialenergy between the particles and potentialenergy between the electrons and the nuclei

    //If we are setting the nuclei side by side we need to make sure that 1/|R| is not included. That is accounted for by other forces
    if(!m_zeroDistance)
        potentialEnergy += charge/solver->trialFunction()->getNucleusDistance();

    //Taking away the electron electron interaction, used for some tests with Hydrogenic wavesfunctions
    if(solver->getElectronInteration())
    {
        // Contribution from electron-electron potential
        double r12 = 0;
        for(int i = 0; i < nParticles; i++) {
            for(int j = i + 1; j < nParticles; j++) {
                r12 = 0;
                for(int k = 0; k < nDimensions; k++) {
                    r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
                }
                potentialEnergy += 1. / sqrt(r12);
            }
        }
    }

    //The interaction between the electrons and the nuclei
    for(int i = 0; i < nParticles; i ++)
    {
            double rip1 = norm(r.row(i) + nucleusHalfDistance.t());
            double rip2 = norm(r.row(i) - nucleusHalfDistance.t());

            potentialEnergy -= charge*( 1./rip1 + 1./rip2  );
    }






    return kineticEnergy + potentialEnergy;
}

double BerylliumTwo::lnDerivativeWaveFunction(const mat &r, VMCSolver *solver)
{
    return 0;
}

double BerylliumTwo::lnSecondDerivativeWaveFunction(const mat &r, VMCSolver *solver)
{
    return 0;
}
