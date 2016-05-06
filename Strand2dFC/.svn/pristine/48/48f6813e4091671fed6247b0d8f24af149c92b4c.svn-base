#include "Strand2dFCSPTurbSA.h"


void Strand2dFCSPTurbSA::lhsVisVariableJacobian(const int& npts,
					     const double* q,
					     const double* qa,
					     double* A)
{
  int iq,iqa,iA;
  double r,u,v,e,a,nu;
  for (int n=0; n<npts; n++){
    iq       = n*nq;
    iqa      = n*nqa;
    iA       = n*nqaGradQa*nq;

    r        = 1./q[iq];
    e        = q[iq+3]*r;
    u        = qa[iqa+1];
    v        = qa[iqa+2];
    nu       = qa[iqa+4];
    a        = r*gm1/rGas;

    A[iA   ] =-r*u;
    A[iA+1 ] = r;
    A[iA+2 ] = 0.;
    A[iA+3 ] = 0.;
    A[iA+4 ] = 0.;

    A[iA+5 ] =-r*v;
    A[iA+6 ] = 0.;
    A[iA+7 ] = r;
    A[iA+8 ] = 0.;
    A[iA+9 ] = 0.;

    A[iA+10] =-a*(e-(u*u+v*v));
    A[iA+11] =-a*u;
    A[iA+12] =-a*v;
    A[iA+13] = a;
    A[iA+14] = 0.;
   
    A[iA+15] =-r*nu;
    A[iA+16] = 0.;
    A[iA+17] = 0.;
    A[iA+18] = 0.;
    A[iA+19] = r;
  }
}
