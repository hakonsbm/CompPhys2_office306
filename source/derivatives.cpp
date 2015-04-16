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
    double ri;

    double derivative = 0;

    for(int i = 0; i < solver->getNParticles(); i++)
    {
        ri = norm(r.row(i));
        derivative += (alpha*alpha - (2.*alpha)/ri);
//        cout << "Particle " << i << " correct : r: " << ri << " derivative: " << (alpha*alpha - (2.*alpha)/ri) << endl;
    }

    return derivative;
}

double Derivatives::analyticalPsi1SDerivative(int particleTag, const mat &r, VMCSolver *solver)
{
    double alpha = solver->getAlpha();
    double r_i = norm(r.row(particleTag));
    double x_i = r(particleTag, 0);
    double y_i = r(particleTag, 1);
    double z_i = r(particleTag, 2);


    double derivative = -alpha*(x_i + y_i + z_i)*exp(-alpha*r_i)/r_i;

    return derivative;
}

double Derivatives::analyticalPsi1SDoubleDerivative(int particleTag, const mat &r, VMCSolver *solver)
{
    double alpha = solver->getAlpha();
    double ri = norm(r.row(particleTag));

    double derivative = alpha * (alpha  - (2./ri) )  * exp(-alpha*ri);

//    cout << "Particle " << particleTag << " New : r: " << ri << " derivative: " << alpha * (alpha - (2./ri) )   << endl;


    return derivative;
}

double Derivatives::analyticalPsi2SDerivative(int particleTag, const mat &r, VMCSolver *solver)
{

    double alpha = solver->getAlpha();
    double r_i = norm(r.row(particleTag));
    double x_i = r(particleTag, 0);
    double y_i = r(particleTag, 1);
    double z_i = r(particleTag, 2);


    double derivative = (1.0L/4.0L)*alpha*(alpha*r_i - 4)*(x_i + y_i + z_i)*exp(-1.0L/2.0L*alpha*r_i)/r_i;

    return derivative;
}

double Derivatives::analyticalPsi2SDoubleDerivative(int particleTag, const mat &r, VMCSolver *solver)
{

//    cout << "ParticleTag : " << particleTag <<endl;
    double alpha = solver->getAlpha();
    double ri = norm(r.row(particleTag));

    double derivative = -1.0L/8.0L*alpha*(pow(alpha, 2)*pow(ri, 2) - 10*alpha*ri + 16)*exp(-1.0L/2.0L*alpha*ri)/ri;

    return derivative;
}

double Derivatives::analyticalPsi2PDerivative(int particleTag, const mat &r, VMCSolver *solver)
{
    double alpha = solver->getAlpha();
    double r_i = norm(r.row(particleTag));
    double x_i = r(particleTag, 0);
    double y_i = r(particleTag, 1);
    double z_i = r(particleTag, 2);

    double derivative = -1.0L/2.0L*alpha*(alpha*pow(x_i, 2) + alpha*x_i*y_i + alpha*x_i*z_i - 2*r_i)*exp(-1.0L/2.0L*alpha*r_i)/r_i;

    return derivative;
}

double Derivatives::analyticalPsi2PDoubleDerivative(int particleTag, int dimension, const mat &r, VMCSolver *solver)
{
    double alpha = solver->getAlpha();
    double r_i = norm(r.row(particleTag)) ;
    double x_i = r(particleTag,dimension);

    double derivative = (1.0L/4.0L)*pow(alpha, 2)*x_i*(alpha*r_i - 8)*exp(-1.0L/2.0L*alpha*r_i)/r_i;

    return derivative;
}


