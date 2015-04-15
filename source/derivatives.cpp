#include "derivatives.h"
#include <vmcsolver.h>

Derivatives::Derivatives()
{

}

Derivatives::~Derivatives()
{

}

double Derivatives::numericalDoubleDerivative(const mat &r, VMCSolver *solver)
{
//    cout << "Doing a numerical derivation" << endl;
    double nParticles = solver->getNParticles();
    double nDimensions = solver->getNDimensions();
    double charge = solver->getCharge();
    double h = solver->getH();
    double h2 = solver->getH2();

    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);

    rPlus = rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    double waveFunctionCurrent = solver->trialFunction()->waveFunction(r, solver);

    //Numerical derivation according to central difference: f''(x)  = [f(x+h) + f(x-h) - 2 f(x)]/h²
    double doubleDerivative = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = solver->trialFunction()->waveFunction(rMinus, solver);
            waveFunctionPlus = solver->trialFunction()->waveFunction(rPlus, solver);
            doubleDerivative -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }
    doubleDerivative = h2 * doubleDerivative;

    return doubleDerivative;
}

double Derivatives::analyticalSimpleDoubleDerivative(const mat &r, VMCSolver *solver)
{
    //This calculates the simple parts of the trialfunctions that are without interaction between the molecules
    //Calculates (nabla Psi_S) /Psi_S
    // d²/dx² (sum_i  e^(-alpha r_i ) ) / sum_i  e^(-alpha r_i )

    double alpha = solver->getAlpha();

    double derivative = 0;

    for(int i = 0; i < solver->getNParticles(); i++)
    {
        derivative += alpha/(norm(r.row(i))) - pow(alpha,2)/2.;
    }

    return derivative;
}

double Derivatives::analyticalPsi1SDerivative(int particleTag, const mat &r, VMCSolver *solver)
{
    double alpha = solver->getAlpha();

    double derivative = - alpha * (r(particleTag,1) + r(particleTag,2) + r(particleTag,3)) /* * exp(-alpha*norm(r.row(particleTag)))*/ / norm(r.row(particleTag));

    return derivative;
}

double Derivatives::analyticalPsi1SDoubleDerivative(int particleTag, const mat &r, VMCSolver *solver)
{
    double alpha = solver->getAlpha();

    double derivative = alpha * (alpha * norm(r.row(particleTag)) - 2) /* * exp(-alpha*norm(r.row(particleTag)))*/ / norm(r.row(particleTag));

    return derivative;
}

double Derivatives::analyticalPsi2SDerivative(int particleTag, const mat &r, VMCSolver *solver)
{
    double alpha = solver->getAlpha();

    double derivative = alpha * (alpha * norm(r.row(particleTag))-3) * (r(particleTag,1) + r(particleTag,2) + r(particleTag,3)) /* * exp(-alpha*norm(r.row(particleTag)))*/ / (2 * norm(r.row(particleTag)));

    return derivative;
}

double Derivatives::analyticalPsi2SDoubleDerivative(int particleTag, const mat &r, VMCSolver *solver)
{

//    cout << "ParticleTag : " << particleTag <<endl;
    double alpha = solver->getAlpha();

    double derivative = - alpha * alpha * (alpha * r(particleTag,0) * r(particleTag,0) + alpha * r(particleTag,1) * r(particleTag,1) + alpha * r(particleTag,2) * r(particleTag,2) + 6 / alpha - 6 * norm(r.row(particleTag))) /* * exp(-alpha*norm(r.row(particleTag)))*/ / (2 * norm(r.row(particleTag)));

    return derivative;
}

double Derivatives::analyticalPsi2PDerivative(int particleTag, const mat &r, VMCSolver *solver)
{
    double alpha = solver->getAlpha();

    double derivative = -alpha * (alpha * norm(r.row(particleTag)) - 2) * (r(particleTag,1) + r(particleTag,2) + r(particleTag,3)) /* * exp(-alpha*norm(r.row(particleTag)) / 2) */ / (2 * norm(r.row(particleTag)));

    return derivative;
}

double Derivatives::analyticalPsi2PDoubleDerivative(int particleTag, const mat &r, VMCSolver *solver)
{
    double alpha = solver->getAlpha();

    double derivative = alpha * alpha * (alpha * r(particleTag,1) * r(particleTag,1) + alpha * r(particleTag,2) * r(particleTag,2) + alpha * r(particleTag,3) * r(particleTag,3) + 8 / alpha - 8 * norm(r.row(particleTag))) /* * exp(-alpha*norm(r.row(particleTag)) / 2)*/ / (4 * norm(r.row(particleTag)));

    return derivative;
}


