#include "StrandSPTurbSA.h"


void StrandSPTurbSA::lhsInvFluxJacobian(const int& npts,
					const double* A,
					const double* xv,
					const double* q,
					const double* qa,
					double* M)
{
  int iq,iqa,iA,iM;
  double Ax,Ay,u,v,qq,qn,h,t,nu;

  for (int n=0; n<npts; n++){
    iq       = nq   *n;
    iqa      = nqa  *n;
    iA       = ndim *n;
    iM       = nq*nq*n;
    Ax       = A[iA  ];
    Ay       = A[iA+1];

    u        = qa[iqa+1];
    v        = qa[iqa+2];
    t        = qa[iqa+3];
    nu       = qa[iqa+4];
    qq       = .5*(u*u+v*v);
    qn       = u*Ax+v*Ay;
    h        = ggm1*rGas*t+qq;

    M[iM   ] =-xv[n];
    M[iM+1 ] = Ax;
    M[iM+2 ] = Ay;
    M[iM+3 ] = 0.;
    M[iM+4 ] = 0.;

    M[iM+5 ] = Ax*gm1*qq-u*qn;
    M[iM+6 ] = u*Ax-gm1*u*Ax+qn-xv[n];
    M[iM+7 ] = u*Ay-gm1*v*Ax;
    M[iM+8 ] = Ax*gm1;
    M[iM+9 ] = 0.;

    M[iM+10] = Ay*gm1*qq-v*qn;
    M[iM+11] = v*Ax-gm1*u*Ay;
    M[iM+12] = v*Ay-gm1*v*Ay+qn-xv[n];
    M[iM+13] = Ay*gm1;
    M[iM+14] = 0.;

    M[iM+15] = qn*(gm1*qq-h);
    M[iM+16] = Ax*h-gm1*u*qn;
    M[iM+17] = Ay*h-gm1*v*qn;
    M[iM+18] = gamma*qn-xv[n];
    M[iM+19] = 0.;
    
    M[iM+20] =-qn*nu;
    M[iM+21] = Ax*nu;
    M[iM+22] = Ay*nu;
    M[iM+23] = 0.;
    M[iM+24] = qn-xv[n];
  }
}
