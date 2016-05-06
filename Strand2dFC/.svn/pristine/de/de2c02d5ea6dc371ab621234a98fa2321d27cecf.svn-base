#include "Strand2dFCSPLam.h"


void Strand2dFCSPLam::rhsVisFluxNCoeff(const int& npts,
				       const double* jac,
				       const double* xs,
				       const double* ys,
				       const double* xn,
				       const double* yn,
				       const double* q,
				       const double* qa,
				       double* b)
{
  int ib,iqa,iqas;
  double u,v,mu,kp,a,b22,b23,b32,b33;
  for (int n=0; n<npts; n++){
    ib       = nq*nqaGradQa*n;
    iqa      = nqa*n;
    iqas     = nqaGradQa*n;
    u        = qa[iqa+1];
    v        = qa[iqa+2];
    mu       = qa[iqa+4];
    kp       = qa[iqa+5];
    a        = xs[n]*xs[n]+ys[n]*ys[n];
    b22      = mu*(a+ys[n]*ys[n]/3.)/jac[n];
    b23      =-mu*(  xs[n]*ys[n]/3.)/jac[n];
    b32      =-mu*(  xs[n]*ys[n]/3.)/jac[n];
    b33      = mu*(a+xs[n]*xs[n]/3.)/jac[n];

    b[ib   ] = 0.;
    b[ib+ 1] = 0.;
    b[ib+ 2] = 0.;

    b[ib+ 3] = b22;
    b[ib+ 4] = b23;
    b[ib+ 5] = 0.;

    b[ib+ 6] = b32;
    b[ib+ 7] = b33;
    b[ib+ 8] = 0.;

    b[ib+ 9] = u*b22+v*b32;
    b[ib+10] = u*b23+v*b33;
    b[ib+11] = a*kp/jac[n];
  }
}
