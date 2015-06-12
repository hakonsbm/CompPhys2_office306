#include "helium.h"

Helium::Helium(VMCSolver *solver)
{

    simpleFlag = false;
    m_analytical = false;

    m_outfileName = "Helium";

    solver->setCharge(2);
    solver->setNParticles(2);
    solver->setAlpha(1.843);
    solver->setBeta(0.34);

    spin << 0 << 1;
}



double Helium::waveFunction(const mat &r, VMCSolver *solver)
{
    double rSingleParticle, SD, rij, a;
    double product = 1.0;
    double alpha = solver -> getAlpha();
    double beta = solver -> getBeta();

    //Calculate the Jastrow factor
    if(solver->getElectronInteration())
    {
        for(int i = 0; i < solver->getNParticles(); i++) {
            rSingleParticle = 0;
            for(int j = i + 1; j < solver->getNParticles(); j++) {
                rij = 0;
                for(int k = 0; k < solver->getNDimensions(); k++) {
                    rij += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
                }
                rij = sqrt(rij);
                a = spinFactor(i,j);
                product = product * exp(a*rij/(1+beta*rij));
            }
        }
    }

//    cout << product << endl;

//    cout << "Before SD" << endl;
    SD = solver->determinant()->calculateDeterminant(r,alpha,solver); //SlaterDeterminant(r, alpha, solver);
//    cout << SD << endl;


    return SD*product;
}

double Helium::localEnergy(const mat &r, VMCSolver *solver)
{

    double kineticEnergy = 0;
    double potentialEnergy = 0;

    double r1 = norm(r.row(0));
    double r2 = norm(r.row(1));
    double r12 = norm(r.row(0) - r.row(1));

    double charge = solver->getCharge();


    if(m_analytical)
    {
    //Calculates the kinetic energy as the ratios of -1/2* ( d²/dx²|D| /|D| + 2 * (d/dx |D|/|D|)*d/dx Psi_C/Psi_C + d²/dx² Psi_C /Psi_C )

        //If we want to compute without electroninteraction with a simplified version ofthe trialfunction containing only the slater determinant
        //then only the slater determinant ratio laplacian is used.
        if(solver->getElectronInteration())
        {
            solver->derivatives()->analyticalLaplacianRatio(kineticEnergy, r, solver);
        }
        else
        {
            kineticEnergy += solver->determinant()->laplacianSlaterDeterminant(r, solver);
        }

            kineticEnergy *= -1./2.;

//            cout << kineticEnergy << " vs " << endl << solver->derivatives()->numericalDoubleDerivative(r, solver) / (2.*waveFunction(r, solver)) << endl;
    }
    else
        kineticEnergy =  solver->derivatives()->numericalDoubleDerivative(r, solver) / (2.*waveFunction(r, solver)); //The minus looks to be baked into the
                                                                                                                     //doubleDerivative Function


    //Taking away the electron electron interaction, used for some tests with Hydrogenic wavesfunctions
    if(solver->getElectronInteration())
        potentialEnergy = - charge*(1./r1+1./r2) + 1./(r12);
    else
        potentialEnergy = - charge*(1./r1+1./r2);

//    cout << "potential is " << -solver->determinant()->laplacianSlaterDeterminant(r, solver)/2. + potentialEnergy << endl;



    return kineticEnergy + potentialEnergy;

}

double Helium::lnDerivativeWaveFunction(const mat &r, VMCSolver *solver)
{
    double r1 = norm(r.row(0));
    double r2 = norm(r.row(1));
    double r12 = norm(r.row(0) - r.row(1));
    double dotproduct = dot(r.row(0),r.row(1));
    double derivative = 0;
    double beta = solver->getBeta();

    derivative = (-1*solver->getAlpha()*(r1+r2)*(1-dotproduct/(r1*r2))*pow(1+r12*beta,2)+3*r12*beta+r12+3)/pow(1+r12*beta,5);

    return derivative;
}

double Helium::lnSecondDerivativeWaveFunction(const mat &r, VMCSolver *solver)
{
    double r1 = norm(r.row(0));
    double r2 = norm(r.row(1));
    double r12 = norm(r.row(0) - r.row(1));
    double dotproduct = dot(r.row(0),r.row(1));
    double factor = solver->getAlpha()*(r1+r2)*(1-dotproduct/(r1*r2));
    double secondDerivative = 0;
    double beta = solver->getBeta();

    secondDerivative = (r12*(3*r12*r12*factor*beta*beta+6*r12*factor*beta-12*r12*beta-5*r12+3*factor-12))/pow(1+r12*beta,6);

    return secondDerivative;
}
