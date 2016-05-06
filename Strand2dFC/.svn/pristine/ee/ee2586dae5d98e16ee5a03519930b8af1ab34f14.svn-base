#include "Strand2dFCSPLam.h"


void Strand2dFCSPLam::rhsVisFluxS(const int& npts,
				  const double* jac,
				  const double* xs,
				  const double* ys,
				  const double* xn,
				  const double* yn,
				  const double* q,
				  const double* qa,
				  const double* qas,
				  double* f)
{
  int iq,iqa,iqas;
  double u,v,mu,kp,us,vs,Ts,a,b22,b23,b32,b33;
  for (int n=0; n<npts; n++){
    iq      = nq *n;
    iqa     = nqa*n;
    iqas    = nqaGradQa*n;
    u       = qa[iqa+1];
    v       = qa[iqa+2];
    mu      = qa[iqa+4];
    kp      = qa[iqa+5];
    us      = qas[iqas  ]; //gradQa numbering
    vs      = qas[iqas+1];
    Ts      = qas[iqas+2];
    a       = xs[n]*xn[n]+ys[n]*yn[n];
    b22     =-mu*( a     +ys[n]*yn[n]/3.)/jac[n];
    b23     = mu*( jac[n]+ys[n]*xn[n]/3.)/jac[n];
    b32     = mu*(-jac[n]+xs[n]*yn[n]/3.)/jac[n];
    b33     =-mu*( a     +xs[n]*xn[n]/3.)/jac[n];
    f[iq  ] = 0.;
    f[iq+1] = b22*us+b23*vs;
    f[iq+2] = b32*us+b33*vs;
    f[iq+3] =(u*b22+v*b32)*us+(u*b23+v*b33)*vs-a*kp*Ts/jac[n];
  }
}
