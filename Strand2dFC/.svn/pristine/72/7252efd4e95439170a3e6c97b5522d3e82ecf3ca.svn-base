#include "Strand2dFCSPLam.h"


void Strand2dFCSPLam::lhsVisVariableJacobian(const int& npts,
					     const double* q,
					     const double* qa,
					     double* A)
{
  int iq,iqa,iA;
  double r,u,v,e,a;
  for (int n=0; n<npts; n++){
    iq       = n*nq;
    iqa      = n*nqa;
    iA       = n*nqaGradQa*nq;

    r        = 1./q[iq];
    e        = q[iq+3]*r;
    u        = qa[iqa+1];
    v        = qa[iqa+2];
    a        = r*gm1/rGas;

    A[iA   ] =-r*u;
    A[iA+1 ] = r;
    A[iA+2 ] = 0.;
    A[iA+3 ] = 0.;

    A[iA+4 ] =-r*v;
    A[iA+5 ] = 0.;
    A[iA+6 ] = r;
    A[iA+7 ] = 0.;

    A[iA+8 ] =-a*(e-(u*u+v*v));
    A[iA+9 ] =-a*u;
    A[iA+10] =-a*v;
    A[iA+11] = a;
  }
}
