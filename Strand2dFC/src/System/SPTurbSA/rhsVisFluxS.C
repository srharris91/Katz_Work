#include "Strand2dFCSPTurbSA.h"
#include <math.h>

void Strand2dFCSPTurbSA::rhsVisFluxS(const int& npts,
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
  double rnu,u,v,nu,mu,k,us,vs,Ts,nus,a,chi,chi3,fv1,mut,kt,fn,qs,nvs,
    b22,b23,b32,b33;

  for (int n=0; n<npts; n++){
    iq      = nq *n;
    iqa     = nqa*n;
    iqas    = nqaGradQa*n;
    rnu     = q  [iq  +4];
    u       = qa [iqa +1];
    v       = qa [iqa +2];
    nu      = qa [iqa +4];
    mu      = qa [iqa +5];
    k       = qa [iqa +6];
    us      = qas[iqas  ]; //gradQa numbering
    vs      = qas[iqas+1];
    Ts      = qas[iqas+2];
    nus     = qas[iqas+3];
    a       = xs[n]*xn[n]+ys[n]*yn[n];    

    chi     = rnu/mu;
    chi3    = chi*chi*chi;
    fv1     = chi3/(chi3+cv1*cv1*cv1);
    mut     = rnu*fv1;
    kt      = mut*rGas*ggm1/PrnT;
    fn      =(mu+rnu)/sigma;
    if (rnu < 0.){
      fn    =(cn1+chi3)/(cn1-chi3);
      fn    =(mu+rnu*fn)/sigma;
      mut   = 0.;
      kt    = 0.;
    }
    mu     += mut;
    k      += kt;

    qs      = k*Ts;
    nvs     = fn*nus; 

    b22     =-mu*( a     +ys[n]*yn[n]/3.)/jac[n];
    b23     = mu*( jac[n]+ys[n]*xn[n]/3.)/jac[n];
    b32     = mu*(-jac[n]+xs[n]*yn[n]/3.)/jac[n];
    b33     =-mu*( a     +xs[n]*xn[n]/3.)/jac[n];

    f[iq  ] = 0.;
    f[iq+1] = b22*us+b23*vs;
    f[iq+2] = b32*us+b33*vs;
    f[iq+3] =(u*b22+v*b32)*us+(u*b23+v*b33)*vs-a*qs/jac[n];;
    f[iq+4] =-a*nvs/jac[n];
  }
}
