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
    double kineticEnergy, potentialEnergy;

    double r1 = norm(r.row(0));
    double r2 = norm(r.row(1));
    double r12 = norm(r.row(0) - r.row(1));
    vec gradientSlater, gradientJastrow;


    double charge = solver->getCharge();


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
        potentialEnergy = - charge*(1./r1+1./r2) + 1./(r12);
    else
        potentialEnergy = - charge*(1./r1+1./r2);



    return kineticEnergy + potentialEnergy;

}
