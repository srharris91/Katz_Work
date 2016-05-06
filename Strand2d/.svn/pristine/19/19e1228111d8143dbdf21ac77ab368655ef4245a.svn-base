#include "StrandSPTurbSA.h"


void StrandSPTurbSA::stepQAdd(const int& npts,
			      const double* q,
			      double* qa)
{
  int iq,iqa;
  int j=1;
  double r,ru,rv,re,u,v,e,p,t,mu,k,nu,dw,rnu;
  for (int n=0; n<npts; n++){
    iq        = nq *n;
    iqa       = nqa*n;
    r         = q[iq  ];
    ru        = q[iq+1];
    rv        = q[iq+2];
    re        = q[iq+3];
    rnu       = q[iq+4];
    u         = ru/r;
    v         = rv/r;
    nu        = rnu/r;
    p         = gm1*(re-.5*r*(u*u+v*v));
    t         = p/(r*rGas);
    transport.getViscosity   (j,&p,&t,&mu);
    transport.getConductivity(j,&p,&t,&k );
    qa[iqa  ] = p;
    qa[iqa+1] = u;
    qa[iqa+2] = v;
    qa[iqa+3] = t;
    qa[iqa+4] = nu;
    qa[iqa+5] = mu;
    qa[iqa+6] = k;
    // wall distance already set
  }
}
