#include "Strand2dFCSPLam.h"


void Strand2dFCSPLam::rhsVisFluxX(const int& npts,
				  const double* q,
				  const double* qa,
				  const double* qax,
				  const double* qay,
				  double* f)
{
  int iq,iqa,iqax;
  double u,v,mu,kp,ux,vx,uy,vy,Tx,dd,sxx,sxy;
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
    dd      =(ux+vy)/3.;
    sxx     = 2.*mu*(ux-dd);
    sxy     =    mu*(uy+vx);
    f[iq  ] = 0.;
    f[iq+1] = sxx;
    f[iq+2] = sxy;
    f[iq+3] = u*sxx+v*sxy+kp*Tx;
  }
}
