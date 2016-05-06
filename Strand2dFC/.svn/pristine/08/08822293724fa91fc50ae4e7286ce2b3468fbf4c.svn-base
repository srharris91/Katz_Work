#include "Strand2dFCSPTurbSA.h"
#include <math.h>

void Strand2dFCSPTurbSA::rhsVisFlux(const int& npts,
			      const double* q,
			      const double* qa,
			      const double* qax,
			      const double* qay,
			      double* f,
			      double* g)
{
  int iq,iqa,iqax;
  double rnu,u,v,nu,mu,k,ux,uy,vx,vy,Tx,Ty,nux,nuy,chi,chi3,fv1,mut,kt,fn,
    dd,sxx,sxy,syx,syy,qxx,qyy,nvx,nvy;

  for (int n=0; n<npts; n++){
    iq      = nq *n;
    iqa     = nqa*n;
    iqax    = nqaGradQa*n;
    rnu     = q  [iq  +4];
    u       = qa [iqa +1];
    v       = qa [iqa +2];
    nu      = qa [iqa +4];
    mu      = qa [iqa +5];
    k       = qa [iqa +6];
    ux      = qax[iqax  ]; //gradQa numbering
    uy      = qay[iqax  ];
    vx      = qax[iqax+1];
    vy      = qay[iqax+1];
    Tx      = qax[iqax+2];
    Ty      = qay[iqax+2];
    nux     = qax[iqax+3];
    nuy     = qay[iqax+3];
 
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

    dd      =(ux+vy)/3.;
    sxx     = 2.*mu*(ux-dd);
    sxy     =    mu*(uy+vx);
    syx     = sxy;
    syy     = 2.*mu*(vy-dd);
    qxx     =-k*Tx;
    qyy     =-k*Ty;
    nvx     = fn*nux;
    nvy     = fn*nuy;

    f[iq  ] = 0.;
    f[iq+1] = sxx;
    f[iq+2] = sxy;
    f[iq+3] = u*sxx+v*sxy-qxx;
    f[iq+4] = nvx;

    g[iq  ] = 0.;
    g[iq+1] = syx;
    g[iq+2] = syy;
    g[iq+3] = u*syx+v*syy-qyy;
    g[iq+4] = nvy;
  }
}
