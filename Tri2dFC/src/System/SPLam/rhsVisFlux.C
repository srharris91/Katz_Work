#include "Tri2dFCSPLam.h"


void Tri2dFCSPLam::rhsVisFlux(const int& npts,
			      const double* q,
			      const double* qa,
			      const double* qax,
			      const double* qay,
			      double* f,
			      double* g)
{
  int iq,iqa,iqax;
  double u,v,mu,kp,ux,vx,uy,vy,Tx,Ty,dd,sxx,sxy,syx,syy;
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
    Tx      = qax[iqax+2];
    Ty      = qay[iqax+2];
    dd      =(ux+vy)/3.;
    sxx     = 2.*mu*(ux-dd);
    sxy     =    mu*(uy+vx);
    syx     = sxy;
    syy     = 2.*mu*(vy-dd);
    f[iq  ] = 0.;
    f[iq+1] = sxx;
    f[iq+2] = sxy;
    f[iq+3] = u*sxx+v*sxy+kp*Tx;
    g[iq  ] = 0.;
    g[iq+1] = syx;
    g[iq+2] = syy;
    g[iq+3] = u*syx+v*syy+kp*Ty;
  }
}
