#include "Strand2dFCSPTurbSA.h"


void Strand2dFCSPTurbSA::rhsVisFluxY(const int& npts,
				  const double* q,
				  const double* qa,
				  const double* qax,
				  const double* qay,
				  double* g)
{
  int iq,iqa,iqax;
  double u,v,mu,kp,ux,vx,uy,vy,Ty,dd,syx,syy;
  for (int n=0; n<npts; n++){
    iq      = nq *n;
    iqa     = nqa*n;
    iqax    = nqaGradQa*n;
    u       = qa[iqa+1];
    v       = qa[iqa+2];
    mu      = qa[iqa+4];
    kp      = qa[iqa+5];
    ux      = qax[iqax  ]; //gradQa numbering
    vx      = qax[iqax+1];
    uy      = qay[iqax  ];
    vy      = qay[iqax+1];
    Ty      = qay[iqax+2];
    dd      =(ux+vy)/3.;
    syx     =    mu*(uy+vx);
    syy     = 2.*mu*(vy-dd);
    g[iq  ] = 0.;
    g[iq+1] = syx;
    g[iq+2] = syy;
    g[iq+3] = u*syx+v*syy+kp*Ty;
  }
}
