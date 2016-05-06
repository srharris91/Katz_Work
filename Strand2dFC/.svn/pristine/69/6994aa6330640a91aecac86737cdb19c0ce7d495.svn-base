#include "Strand2dFCSPTurbSA.h"
#include <math.h>

void Strand2dFCSPTurbSA::rhsVisFluxNCoeff(const int& npts,
				       const double* jac,
				       const double* xs,
				       const double* ys,
				       const double* xn,
				       const double* yn,
				       const double* q,
				       const double* qa,
				       double* b)
{
  int ib,iq,iqa;
  double rnu,u,v,nu,mu,k,a,chi,chi3,fv1,mut,kt,fn,b22,b23,b32,b33;
  for (int n=0; n<npts; n++){
    ib       = nq*nqaGradQa*n;
    iq       = nq*n;
    iqa      = nqa*n;
    rnu      = q [iq +4];
    u        = qa[iqa+1];
    v        = qa[iqa+2];
    nu       = qa[iqa+4];
    mu       = qa[iqa+5];
    k        = qa[iqa+6];    
    a        = xs[n]*xs[n]+ys[n]*ys[n];
  
    chi      = rnu/mu;
    chi3     = chi*chi*chi;
    fv1      = chi3/(chi3+cv1*cv1*cv1);
    mut      = rnu*fv1;
    kt       = mut*rGas*ggm1/PrnT;
    fn       =(mu+rnu)/sigma;
    if (rnu < 0.){
      fn     =(cn1+chi3)/(cn1-chi3);
      fn     =(mu+rnu*fn)/sigma;
      mut    = 0.;
      kt     = 0.;
    }
    mu      += mut;
    k       += kt;

    b22      = mu*(a+ys[n]*ys[n]/3.)/jac[n];
    b23      =-mu*(  xs[n]*ys[n]/3.)/jac[n];
    b32      =-mu*(  xs[n]*ys[n]/3.)/jac[n];
    b33      = mu*(a+xs[n]*xs[n]/3.)/jac[n];

    b[ib   ] = 0.;
    b[ib+ 1] = 0.;
    b[ib+ 2] = 0.;
    b[ib+ 3] = 0.;

    b[ib+ 4] = b22;
    b[ib+ 5] = b23;
    b[ib+ 6] = 0.;
    b[ib+ 7] = 0.;

    b[ib+ 8] = b32;
    b[ib+ 9] = b33;
    b[ib+10] = 0.;
    b[ib+11] = 0.;

    b[ib+12] = u*b22+v*b32;
    b[ib+13] = u*b23+v*b33;
    b[ib+14] = a*k/jac[n];
    b[ib+15] = 0.;

    b[ib+16] = 0.;
    b[ib+17] = 0.;
    b[ib+18] = 0.;
    b[ib+19] = a*fn/jac[n];
  }
}
