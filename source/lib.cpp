#include "lib.h"
#include <cmath>
#include <armadillo>
#include <new>
using namespace arma;
/*
** The function
**         ran2()
** is a long periode (> 2 x 10^18) random number generator of
** L'Ecuyer and Bays-Durham shuffle and added safeguards.
** Call with idum a negative integer to initialize; thereafter,
** do not alter idum between sucessive deviates in a
** sequence. RNMX should approximate the largest floating point value
** that is less than 1.
** The function returns a uniform deviate between 0.0 and 1.0
** (exclusive of end-point values).
*/

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
  int            j;
  long           k;
  static long    idum2 = 123456789;
  static long    iy=0;
  static long    iv[NTAB];
  double         temp;

  if(*idum <= 0) {
    if(-(*idum) < 1) *idum = 1;
    else             *idum = -(*idum);
    idum2 = (*idum);
    for(j = NTAB + 7; j >= 0; j--) {
      k     = (*idum)/IQ1;
      *idum = IA1*(*idum - k*IQ1) - k*IR1;
      if(*idum < 0) *idum +=  IM1;
      if(j < NTAB)  iv[j]  = *idum;
    }
    iy=iv[0];
  }
  k     = (*idum)/IQ1;
  *idum = IA1*(*idum - k*IQ1) - k*IR1;
  if(*idum < 0) *idum += IM1;
  k     = idum2/IQ2;
  idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
  if(idum2 < 0) idum2 += IM2;
  j     = iy/NDIV;
  iy    = iv[j] - idum2;
  iv[j] = *idum;
  if(iy < 1) iy += IMM1;
  if((temp = AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

// End: function ran2()


// random numbers with gaussian distribution
double GaussianDeviate(long * idum)
{
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;

  if ( idum < 0) iset =0;
  if (iset == 0) {
    do {
      v1 = 2.*ran2(idum) -1.0;
      v2 = 2.*ran2(idum) -1.0;
      rsq = v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.);
    fac = sqrt(-2.*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  } else {
    iset =0;
    return gset;
  }
} // end function for gaussian deviates

// the function ludcmp from c++ library file from FYS3150, modified to use armadillo
void ludcmp(mat &a, int n, int *indx, double *d)
{
   int      i, imax, j, k;
   double   big, dum, sum, temp, *vv;

  vv = new(std::nothrow) double [n];
  if(!vv) {
    printf("\n\nError in function ludcm():");
    printf("\nNot enough memory for vv[%d]\n",n);
    exit(1);
  }

   *d = 1.0;                              // no row interchange yet
   for(i = 0; i < n; i++) {     // loop over rows to get scaling information
      big = ZERO;
      for(j = 0; j < n; j++) {
         if((temp = fabs(a(i, j))) > big) big = temp;
      }
      if(big == ZERO) {
         printf("\n\nSingular matrix in routine ludcmp()\n");
         exit(1);
      }
      vv[i] = 1.0/big;                 // save scaling */
   } // end i-loop */

   for(j = 0; j < n; j++) {     // loop over columns of Crout's method
      for(i = 0; i< j; i++) {   // not i = j
         sum = a(i, j);
     for(k = 0; k < i; k++) sum -= a(i, k)*a(k, j);
     a(i, j) = sum;
      }
      big = ZERO;   // initialization for search for largest pivot element
      for(i = j; i< n; i++) {
         sum = a(i, j);
     for(k = 0; k < j; k++) sum -= a(i, k)*a(k, j);
     a(i, j) = sum;
     if((dum = vv[i]*fabs(sum)) >= big) {
        big = dum;
        imax = i;
     }
      } // end i-loop
      if(j != imax) {    // do we need to interchange rows ?
         for(k = 0;k< n; k++) {       // yes
        dum        = a(imax, k);
        a(imax, k) = a(j, k);
        a(j, k)    = dum;
     }
     (*d)    *= -1;            // and change the parit of d
     vv[imax] = vv[j];         // also interchange scaling factor
      }
      indx[j] = imax;
      if(fabs(a(j, j)) < ZERO)  a(j, j) = ZERO;

        /*
        ** if the pivot element is zero the matrix is singular
        ** (at least to the precision of the algorithm). For
        ** some application of singular matrices, it is desirable
        ** to substitute ZERO for zero,
        */

      if(j < (n - 1)) {                   // divide by pivot element
         dum = 1.0/a(j, j);
     for(i=j+1;i < n; i++) a(i, j) *= dum;
      }
   } // end j-loop over columns

   delete [] vv;   // release local memory

}  // End: function ludcmp()
