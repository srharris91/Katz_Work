#include "StrandSPLam.h"


void StrandSPLam::lhsInvFluxJacobian(const int& npts,
				     const double* A,
				     const double* xv,
				     const double* q,
				     const double* qa,
				     double* M)
{
  int iq,iqa,iA,iM;
  double Ax,Ay,u,v,qq,qn,h;

  for (int n=0; n<npts; n++){
    iq      = nq   *n;
    iqa     = nqa  *n;
    iA      = ndim *n;
    iM      = nq*nq*n;
    Ax      = A[iA  ];
    Ay      = A[iA+1];

    u       = qa[iqa+1];
    v       = qa[iqa+2];
    qq      = .5*(u*u+v*v);
    qn      = u*Ax+v*Ay;
    h       = ggm1*rGas*qa[iqa+3]+qq;

    M[iM   ] =-xv[n];
    M[iM+1 ] = Ax;
    M[iM+2 ] = Ay;
    M[iM+3 ] = 0.;

    M[iM+4 ] = Ax*gm1*qq-u*qn;
    M[iM+5 ] = u*Ax-gm1*u*Ax+qn-xv[n];
    M[iM+6 ] = u*Ay-gm1*v*Ax;
    M[iM+7 ] = Ax*gm1;

    M[iM+8 ] = Ay*gm1*qq-v*qn;
    M[iM+9 ] = v*Ax-gm1*u*Ay;
    M[iM+10] = v*Ay-gm1*v*Ay+qn-xv[n];
    M[iM+11] = Ay*gm1;

    M[iM+12] = qn*(gm1*qq-h);
    M[iM+13] = Ax*h-gm1*u*qn;
    M[iM+14] = Ay*h-gm1*v*qn;
    M[iM+15] = gamma*qn-xv[n];
  }
}
