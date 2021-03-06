#include "StrandSPLam.h"


void StrandSPLam::stepQAdd(const int& npts,
			   const double* q,
			   double* qa)
{
  int iq,iqa;
  int j=1;
  double r,ru,rv,re,u,v,e,p,t,mu,k;
  for (int n=0; n<npts; n++){
    iq        = nq*n;
    iqa       = nqa*n;
    r         = q[iq  ];
    ru        = q[iq+1];
    rv        = q[iq+2];
    re        = q[iq+3];
    u         = ru/r;
    v         = rv/r;
    p         = gm1*(re-.5*r*(u*u+v*v));
    t         = p/(r*rGas);
    transport.getViscosity   (j,&p,&t,&mu);
    transport.getConductivity(j,&p,&t,&k );
    qa[iqa  ] = p;
    qa[iqa+1] = u;
    qa[iqa+2] = v;
    qa[iqa+3] = t;
    qa[iqa+4] = mu;
    qa[iqa+5] = k;
  }
}
