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

vec Derivatives::analyticalPsi1SDerivative(int particleTag, const mat &r, VMCSolver *solver)
{
    double alpha = solver->getAlpha();
    double r_i = norm(r.row(particleTag));


    vec derivative = -alpha*r.row(particleTag)*exp(-alpha*r_i)/r_i;

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

vec Derivatives::analyticalPsi2SDerivative(int particleTag, const mat &r, VMCSolver *solver)
{

    double alpha = solver->getAlpha();
    double r_i = norm(r.row(particleTag));


    vec derivative = (1.0L/4.0L)*alpha*(alpha*r_i - 4)*r.row(particleTag)*exp(-1.0L/2.0L*alpha*r_i)/r_i;

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

vec Derivatives::analyticalPsi2PDerivative(int particleTag, int dimension, const mat &r, VMCSolver *solver)
{
    double alpha = solver->getAlpha();
    double r_i = norm(r.row(particleTag));
    double x_i = r(particleTag, 0);
    double y_i = r(particleTag, 1);
    double z_i = r(particleTag, 2);
    double o_i = r(particleTag, dimension);
    double factor = -1.0L/2.0L*o_i*alpha*exp(-1.0L/2.0L*alpha*r_i)/r_i;
    vec derivative;

    derivative(0) = factor*(alpha*x_i - 2*r_i);
    derivative(1) = factor*(alpha*y_i - 2*r_i);
    derivative(2) = factor*(alpha*z_i - 2*r_i);

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

double Derivatives::fDerivative(int i, int j, const mat &r, VMCSolver *solver)
{
    //Calculates the d/dx f_ij derivative
    double beta = solver->getBeta();
    double a = solver->trialFunction()->spinFactor(i,j);
    double rij = norm(r.row(i) - r.row(j));

    return a/pow(1+beta*rij, 2);

}

double Derivatives::fDoubleDerivative(int i, int j, const mat &r, VMCSolver *solver)
{
    //Calculates the d²/dx² f_ij derivative
    double beta = solver->getBeta();
    double a = solver->trialFunction()->spinFactor(i,j);
    double rij = norm(r.row(i) - r.row(j));

    return -2*a*beta/pow(1+beta*rij, 3);
}



vec Derivatives::analyticalCorrelationDerivative( const mat &r, VMCSolver *solver)
{
    //This sums over all the electrons and calculates the total correlation gradient ratio term

    int nParticles = solver->getNParticles();
    vec gradient = zeros (3);
    vec rik = zeros (3);
    vec rki = zeros (3);

    //Calculates the interaction from all the particles earlier
    for(int k = 0; k < nParticles; k++)
    {
           for(int i = 0; i < k; i ++)
           {
               rik = (r.row(i) - r.row(k)).t();
               gradient = gradient + rik / norm(rik) * fDerivative(i,k,r, solver);

           }
           for(int i = k +1 ; i < nParticles  ; i ++)
           {
               rki = (r.row(k) - r.row(i)).t();
               gradient -= rki / norm(rki) * fDerivative(k,i,r,solver);
               //cout << norm(rki) << endl;
           }
    }


    return gradient;
}

double Derivatives::analyticalCorrelationDoubleDerivative(const mat &r, VMCSolver *solver)
{
    //Not properly tested yet

    vec rki = zeros(3);
    vec rkj = zeros(3);

    int nParticles = solver->getNParticles();

    int i, j, k;

//    cout << "Got ehre" << endl;
    double laplacian = 0;

    //Summing over all the electrons
    for (k = 0; k < nParticles ; k ++)
    {
        //Summing over the part with correlation between the other electrons than electron k
        for(i = k + 1 ; i < nParticles ; i ++)
        {
            for(j = i + 1; j < nParticles ; j ++)
            {
//                cout << "Sum 1" << endl;
                rkj = (r.row(k) - r.row(j)).t();
                rki = (r.row(k) - r.row(i)).t();
                laplacian = laplacian + dot(rkj,rkj) / (norm(rki)*norm(rkj)) * fDerivative(k,j,r,solver) * fDerivative(k,i,r,solver) ;

            }
        }

        //Summing over the d²/dx² part of the expression
        for(j = 0; j < nParticles ; j ++)
        {
//            if(j != k)
//            {
//            cout << "Sum 2" << endl;
//            rkj = (r.row(k) - r.row(j)).t();
//            laplacian = laplacian + fDoubleDerivative(k,j, r, solver) + 2.*fDerivative(k,j,r, solver)/ norm(rkj);

//            cout << "Sum 2 end" << endl;
//            }
        }

    }



    return laplacian;
}

//double Derivatives::analyticalCorrelationDoubleDerivative( const mat &r, VMCSolver *solver)
//{
//    //This sums over all the electrons and calculates the total correlation gradient ratio term

//    int nParticles = solver->getNParticles();
//    int nDimensions = solver->getNDimensions();
//    double gradient = 0;
//    vec rik = zeros (3);
//    vec rki = zeros (3);

//    //Calculates the interaction from all the particles earlier
//    for(int k = 0; k < nParticles; k++)
//    {
//        cout << k << endl;
//           for(int i = 0; i < k-1; i ++)
//           {
//               rik = (r.row(i) - r.row(k)).t();
//               gradient += (nDimensions - 1) / norm(rik) * fDerivative(i,k,r,solver) + fDoubleDerivative(i,k,r,solver);
//           }
//           for(int i = k +1 ; i < nParticles - 1 ; i ++)
//           {

//               rki = (r.row(k) - r.row(i)).t();
//               gradient -= (nDimensions - 1) / norm(rki) * fDerivative(k,i,r,solver) + fDoubleDerivative(k,i,r,solver);
//           }
//    }
//    return gradient;//DON'T FORGET TO ADD THE SUMMAND OF (∇Ψ_C/Ψ_C)² TO THE VARIABLE GRADIENT HERE WHEN YOU CALL THE FUNCTION TO GET ∇Ψ/Ψ, SEE EQUATION (16.38)
//}
