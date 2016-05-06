#include "StrandSPLam.h"


void StrandSPLam::lhsPrimVarJacobian(const int& npts,
				     const double* q,
				     const double* qa,
				     double* M)
{
  int iq,iqa,iM;
  double rr,e,u,v,qq,grr;

  for (int n=0; n<npts; n++){
    iq   = nq*n;
    iqa  = nqa*n;
    iM   = nq*nq*n;

    rr    = 1./q[iq];
    e     = q[iq+3]*rr;
    u     = qa[iqa+1];
    v     = qa[iqa+2];
    qq    = u*u+v*v;
    grr   = gm1*rr/rGas;

    M[iM   ] = gm1*qq*.5;
    M[iM+1 ] =-gm1*u;
    M[iM+2 ] =-gm1*v;
    M[iM+3 ] = gm1;
    
    M[iM+4 ] =-rr*u;
    M[iM+5 ] = rr;
    M[iM+6 ] = 0.;
    M[iM+7 ] = 0.;
    
    M[iM+8 ] =-rr*v;
    M[iM+9 ] = 0.;
    M[iM+10] = rr;
    M[iM+11] = 0.;
    
    M[iM+12] = grr*(qq-e);
    M[iM+13] =-grr*u;
    M[iM+14] =-grr*v;
    M[iM+15] = grr;
  }
}
