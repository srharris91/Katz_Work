#include "StrandSPTurbSA.h"
#include <math.h>


void StrandSPTurbSA::lhsSourceJacobianCoarse(const int& npts,
					     const double* v,
					     const double* q,
					     const double* qa,
					     double* M)
{
  int iM;
  for (int n=0; n<npts; n++){
    iM       = nq*nq*n;

    M[iM   ] = 0.;
    M[iM+1 ] = 0.;
    M[iM+2 ] = 0.;
    M[iM+3 ] = 0.;
    M[iM+4 ] = 0.;

    M[iM+5 ] = 0.;
    M[iM+6 ] = 0.;
    M[iM+7 ] = 0.;
    M[iM+8 ] = 0.;
    M[iM+9 ] = 0.;

    M[iM+10] = 0.;
    M[iM+11] = 0.;
    M[iM+12] = 0.;
    M[iM+13] = 0.;
    M[iM+14] = 0.;

    M[iM+15] = 0.;
    M[iM+16] = 0.;
    M[iM+17] = 0.;
    M[iM+18] = 0.;
    M[iM+19] = 0.;

    M[iM+20] = 0.;
    M[iM+21] = 0.;
    M[iM+22] = 0.;
    M[iM+23] = 0.;
    M[iM+24] = 0.;
  }
}
