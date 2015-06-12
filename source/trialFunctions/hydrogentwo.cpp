#include "hydrogentwo.h"

HydrogenTwo::HydrogenTwo(VMCSolver *solver)
{
    simpleFlag = false;
    m_analytical = false;

    m_outfileName = "HydrogenTwo";

//    cout << "got here" << endl;


//    solver->trialFunction()->setNucleusDistance(1.4);

//    cout << "got here" << endl;


    solver->setCharge(2);
    solver->setNParticles(2);
    solver->setAlpha(1.289);
    solver->setBeta(0.401);
//    solver->trialFunction()->setNucleusDistance(1.4);

    spin << 0 << 1;

}

HydrogenTwo::~HydrogenTwo()
{

}

void HydrogenTwo::calculateAlpha(VMCSolver *solver)
{
//    solver->setAlpha(4.);


    return;
}

double HydrogenTwo::waveFunction(const mat &r, VMCSolver *solver)
{
    double rip1, rip2 , SD, rij, a;     //rip1 is the distance from electron i to nucleus 1
    double product = 1.0;
    double alpha = solver -> getAlpha();
    double beta = solver -> getBeta();
    double nParticles = solver->getNParticles();
    double nDimensions = solver->getNDimensions();
    double RHalf = solver->trialFunction()->getNucleusDistance()/2.;
    double singleParticleWavefunction = 1.0;

    //Calculate the single particle functions
    //We'll assume that the R vector is aligned along the z-axis. So the nucleuses are seperated on the z-axis
    for(int i = 0; i < nParticles ; i++)
    {
        rip1 = 0;
        rip2 = 0;

        rip1 = r(i,0)*r(i,0) + r(i,1)*r(i,1) + (r(i,2) + RHalf)*(r(i,2) + RHalf);
        rip1 = sqrt(rip1);

        rip2 = r(i,0)*r(i,0) + r(i,1)*r(i,1) + (r(i,2) - RHalf)*(r(i,2) - RHalf);
        rip2 = sqrt(rip2);

        singleParticleWavefunction *= (exp(-alpha * rip1) + exp(-alpha * rip2) );

        //Pretty sure it is supposed to be divided by two as well, if we want to reproduce the same numbers as Helium when R = 0;
        singleParticleWavefunction /= 2.;

    }

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

//    cout << "Before SD" << endl;
//    SD = solver->determinant()->calculateDeterminant(r,alpha,solver); //SlaterDeterminant(r, alpha, solver);
//    cout << SD << endl;

    return singleParticleWavefunction*product;
}

double HydrogenTwo::localEnergy(const mat &r, VMCSolver *solver)
{
    double kineticEnergy, potentialEnergy;
    int nDimensions = solver->getNDimensions();

    vec nucleusHalfDistance = zeros(nDimensions);
    nucleusHalfDistance(nDimensions - 1) = solver->trialFunction()->getNucleusDistance()/2.;


    double r1p1 = norm(r.row(0) + nucleusHalfDistance.t());
    double r1p2 = norm(r.row(0) - nucleusHalfDistance.t());
    double r2p1 = norm(r.row(1) + nucleusHalfDistance.t());
    double r2p2 = norm(r.row(1) - nucleusHalfDistance.t());


    double r12 = norm(r.row(0) - r.row(1));
    vec gradientSlater, gradientJastrow;


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

//    cout << kineticEnergy << endl;

    //Taking away the electron electron interaction, used for some tests with Hydrogenic wavesfunctions
    if(solver->getElectronInteration())
        potentialEnergy = - (1./r1p1+1./r1p2 + 1./r2p1 + 1./r2p2) + 1./(r12) /*+ 1./norm(2*nucleusHalfDistance)*/;
    else
        potentialEnergy = - (1./r1p1+1./r1p2 + 1./r2p1+1./r2p2);

    if(!m_zeroDistance)
        potentialEnergy += 1./norm(2*nucleusHalfDistance);




    return kineticEnergy + potentialEnergy;
}

double HydrogenTwo::lnDerivativeWaveFunction(const mat &r, VMCSolver *solver)
{
    return 0;
}

double HydrogenTwo::lnSecondDerivativeWaveFunction(const mat &r, VMCSolver *solver)
{
    return 0;
}


