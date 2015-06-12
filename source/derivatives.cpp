#include "derivatives.h"
#include <vmcsolver.h>

Derivatives::Derivatives()
{

}

Derivatives::~Derivatives()
{

}

void Derivatives::numericalGradient(mat &gradient, const mat &r, VMCSolver *solver)
{
    //This calculates a numerical derivative, which is a half of the quantum force

    double nParticles = solver->getNParticles();
    double nDimensions = solver->getNDimensions();
    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);
    double h = solver->getH();

    rPlus = rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    double waveFunctionCurrent = solver->trialFunction()->waveFunction(r, solver);

    // Kinetic energy

//    double kineticEnergy = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = solver->trialFunction()->waveFunction(rMinus, solver);
            waveFunctionPlus = solver->trialFunction()->waveFunction(rPlus, solver);
            gradient(i,j) =  (waveFunctionPlus-waveFunctionMinus);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }

    gradient *= 1/(2.);

}


double Derivatives::numericalDoubleDerivative(const mat &r, VMCSolver *solver)
{
//    cout << "Doing a numerical derivation" << endl;
    double nParticles = solver->getNParticles();
    double nDimensions = solver->getNDimensions();
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

            solver->determinant()->updateSlaterMatrices(rMinus,solver);     // This needs to be updated and changed back at the end to make sure that D(x_old) is correct
            waveFunctionMinus = solver->trialFunction()->waveFunction(rMinus, solver);

            solver->determinant()->updateSlaterMatrices(rPlus,solver);
            waveFunctionPlus = solver->trialFunction()->waveFunction(rPlus, solver);

            doubleDerivative -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }
    solver->determinant()->updateSlaterMatrices(r,solver);   //Returns the slater determinant to the correct position, according to where the particle is now
//    cout << doubleDerivative << endl;
    doubleDerivative = h2 * doubleDerivative;

    return doubleDerivative;
}


void Derivatives::analyticalGradient(mat &gradient, const mat &r, VMCSolver *solver)
{

    gradient += solver->determinant()->gradientSlaterDeterminant(r, solver);

    gradient += solver->derivatives()->analyticalCorrelationGradient( r ,solver);

    //Testing
    //    cout << analyticalPsi1SDerivative(1,r,solver) << endl;
}

void Derivatives::analyticalLaplacianRatio(double &laplacianRatio, const mat &r, VMCSolver *solver)
{
    //This calculates the laplacian ratio of the trialFunction
    // ( d²/dx²|D| /|D| + 2 * (d/dx |D|/|D|)*d/dx Psi_C/Psi_C + d²/dx² Psi_C /Psi_C )

    double tempTerm = 0;
    int nParticles = solver->getNParticles();
    int nDimensions = solver->getNDimensions();
    mat correlationGradient = zeros (nParticles, nDimensions);
    mat slaterGradient = zeros(nParticles, nDimensions);

    laplacianRatio += solver->determinant()->laplacianSlaterDeterminant(r, solver); // d²/dx² Psi_C /Psi_C
//    laplacianRatio += solver->derivatives()->analyticalCorrelationDoubleDerivative(r,solver); //d²/dx²|D| /|D|
    laplacianRatio += analyticalCorrelationLaplacian(r, solver);

    correlationGradient = solver->derivatives()->analyticalCorrelationGradient(r, solver);
    slaterGradient = solver->determinant()->gradientSlaterDeterminant(r, solver);



    //Computes the dot product over all the particles
    for(int i = 0; i < nParticles;i ++)
    {
        tempTerm += dot(slaterGradient.row(i), correlationGradient.row(i) );
    }

    laplacianRatio += 2.*tempTerm; //( d/dx |D|/|D| ) dot d/dx Psi_C/Psi_C


//    cout << "Am in LaplcacianCalculation" << endl;

//////    cout << "double slater is " << solver->determinant()->laplacianSlaterDeterminant(r, solver) << endl;
//    cout << "The correlation d²/dx² is  " << solver->derivatives()->analyticalCorrelationDoubleDerivative(r,solver) << endl;
//    cout << "New attempt at correlation part is: " << analyticalCorrelationLaplacian(r, solver) << endl;
//    cout << "The combined part is  " << 2.*tempTerm << endl;

//    cout << "correlationGradient"  << endl << correlationGradient << endl;
//    cout << "slaterGradient"  << endl << slaterGradient << endl;

//    cout << "TotalEnergy2 is " << -(solver->derivatives()->analyticalCorrelationDoubleDerivative(r,solver) + 2.*tempTerm)/2. << endl;
//    cout << -laplacianRatio/2. << endl;

    return;
}

vec Derivatives::analyticalPsi1SDerivative(int particleTag, const mat &r, VMCSolver *solver)
{

    double alpha = solver->getAlpha();
    double r_i = norm(r.row(particleTag));
    vec derivative = zeros (solver->getNDimensions());

    derivative = (-alpha*r.row(particleTag)*exp(-alpha*r_i)/r_i).t();


    return derivative;
}

double Derivatives::analyticalPsi1SDoubleDerivative(int particleTag, const mat &r, VMCSolver *solver)
{

    double alpha = solver->getAlpha();
    double ri = norm(r.row(particleTag));

    double derivative = alpha * (alpha  - (2./ri) )  * exp(-alpha*ri);



    return derivative;
}

vec Derivatives::analyticalPsi2SDerivative(int particleTag, const mat &r, VMCSolver *solver)
{

    double alpha = solver->getAlpha();
    double r_i = norm(r.row(particleTag));


    vec derivative = ((1.0L/4.0L)*alpha*(alpha*r_i - 4)*r.row(particleTag)*exp(-1.0L/2.0L*alpha*r_i)/r_i).t();

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
    double factor = exp(-0.5*alpha*r_i)/r_i;
    double x_i = r(particleTag,dimension);
    vec gradient = zeros( solver->getNDimensions());
    vec unitVector = zeros (solver->getNDimensions());

    unitVector(dimension) = 1.;

    gradient = factor * (unitVector - 0.5*alpha*x_i*r(particleTag)/r_i);


    return gradient;
}

double Derivatives::analyticalPsi2PDoubleDerivative(int particleTag, int dimension, const mat &r, VMCSolver *solver)
{
    double alpha = solver->getAlpha();
    double r_i = norm(r.row(particleTag)) ;
    double x_i = r(particleTag,dimension);

    double derivative = 0.25*alpha/**alpha*/*x_i*(alpha*r_i - 8)*exp(-0.5*alpha*r_i)/r_i;

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



mat Derivatives::analyticalCorrelationGradient( const mat &r, VMCSolver *solver)
{
    //This sums over all the electrons and calculates the total correlation gradient ratio term

    int nParticles = solver->getNParticles();
    int nDimensions = solver->getNDimensions();
    mat gradient = zeros(nParticles, nDimensions);
    vec rik = zeros (nDimensions);
    vec rki = zeros (nDimensions);


    //Calculates the interaction from all the particles earlier
    for(int k = 0; k < nParticles; k++)
    {
           for(int i = 0; i < k; i ++)
           {
               rik = (r.row(k) - r.row(i)).t();
               gradient.row(k) += (rik / norm(rik) * fDerivative(i,k,r, solver)).t();
           }
           for(int i = k + 1 ; i < nParticles  ; i ++)
           {
               rki = (r.row(i) - r.row(k)).t();
               gradient.row(k) -= (rki / norm(rki) * fDerivative(k,i,r,solver)).t();
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
                rkj = (r.row(j) - r.row(k)).t();
                rki = (r.row(i) - r.row(k)).t();
                laplacian = laplacian + dot(rki,rkj) / (norm(rki)*norm(rkj)) * fDerivative(k,j,r,solver) * fDerivative(k,i,r,solver) ;
            }
        }

        //Summing over the d²/dx² part of the expression
        for(j = 0; j < nParticles ; j ++)
        {
            if(j != k)
            {
            rkj = (r.row(j) - r.row(k)).t();
            laplacian +=  2.*fDerivative(k,j,r, solver)/ norm(rkj) + fDoubleDerivative(k,j, r, solver);
            }
        }

    }

    return laplacian;
}

double Derivatives::analyticalCorrelationLaplacian( const mat &r, VMCSolver *solver)
{
    //This sums over all the electrons and calculates the total correlation gradient ratio term

    int nParticles = solver->getNParticles();
    int nDimensions = solver->getNDimensions();
    double laplacianRatio = 0;
    vec rik = zeros (nDimensions);
    vec rki = zeros (nDimensions);
    mat gradient = zeros(nParticles, nDimensions);


    //Adding in the square of the gradient-correlation-ratio
    gradient = analyticalCorrelationGradient(r,solver);

    //Calculates the interaction from all the particles earlier
    for(int k = 0; k < nParticles; k++)
    {
        laplacianRatio += dot(gradient.row(k),gradient.row(k));
           for(int i = 0; i < k; i ++)
           {
               rik = (r.row(i) - r.row(k)).t();
               laplacianRatio += (nDimensions - 1) / norm(rik) * fDerivative(i,k,r,solver) + fDoubleDerivative(i,k,r,solver);
           }
           for(int i = k + 1 ; i < nParticles ; i ++)
           {
               rki = (r.row(k) - r.row(i)).t();
               laplacianRatio += (nDimensions - 1) / norm(rki) * fDerivative(k,i,r,solver) + fDoubleDerivative(k,i,r,solver);
           }
    }
//    cout << test << endl;
    return laplacianRatio;//DON'T FORGET TO ADD THE SUMMAND OF (∇Ψ_C/Ψ_C)² TO THE VARIABLE GRADIENT HERE WHEN YOU CALL THE FUNCTION TO GET ∇Ψ/Ψ, SEE EQUATION (16.38)
}
