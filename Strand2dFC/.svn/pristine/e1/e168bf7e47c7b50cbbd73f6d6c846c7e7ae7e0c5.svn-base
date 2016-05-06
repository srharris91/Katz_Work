#include "Strand2dFCSPTurbSA.h"


void Strand2dFCSPTurbSA::lhsInvFluxJacobian(const int& npts,
				       const double* Ax,
				       const double* Ay,
				       const double* q,
				       const double* qa,
				       double* A)
{
  int iq,iqa,iA;
  double u,v,qq,qn,h,nu,xv=0.;
  for (int n=0; n<npts; n++){
    iq      = n*nq;
    iqa     = n*nqa;
    iA      = n*nq*nq;

    u       = qa[iqa+1];
    v       = qa[iqa+2];
    nu      = qa[iqa+4];
    qq      = .5*(u*u+v*v);
    qn      = u*Ax[n]+v*Ay[n];
    h       = ggm1*rGas*qa[iqa+3]+qq;

    A[iA   ] =-xv;
    A[iA+1 ] = Ax[n];
    A[iA+2 ] = Ay[n];
    A[iA+3 ] = 0.;
    A[iA+4 ] = 0.;

    A[iA+5 ] = Ax[n]*gm1*qq-u*qn;
    A[iA+6 ] = u*Ax[n]-gm1*u*Ax[n]+qn-xv;
    A[iA+7 ] = u*Ay[n]-gm1*v*Ax[n];
    A[iA+8 ] = Ax[n]*gm1;
    A[iA+9 ] = 0.;

    A[iA+10] = Ay[n]*gm1*qq-v*qn;
    A[iA+11] = v*Ax[n]-gm1*u*Ay[n];
    A[iA+12] = v*Ay[n]-gm1*v*Ay[n]+qn-xv;
    A[iA+13] = Ay[n]*gm1;
    A[iA+14] = 0.;

    A[iA+15] = qn*(gm1*qq-h);
    A[iA+16] = Ax[n]*h-gm1*u*qn;
    A[iA+17] = Ay[n]*h-gm1*v*qn;
    A[iA+18] = gamma*qn-xv;
    A[iA+19] = 0.;
   
    A[iA+20] =-qn*nu;
    A[iA+21] = Ax[n]*nu;
    A[iA+22] = Ay[n]*nu;
    A[iA+23] = 0.;
    A[iA+24] = qn-xv;

  }
}
