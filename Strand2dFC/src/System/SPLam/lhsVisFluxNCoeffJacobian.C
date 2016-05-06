#include "Strand2dFCSPLam.h"


void Strand2dFCSPLam::lhsVisFluxNCoeffJacobian(const int& npts,
					       const double* jac,
					       const double* xs,
					       const double* ys,
					       const double* xn,
					       const double* yn,
					       const double* q,
					       const double* qa,
					       const double* qi,
					       const double* qai,
					       double* A)
{
  int iq,iqa,iA;
  double r,e,u,v,mu,ui,vi,Ti,a,b,c,d,p,t,f,dmudT,dkdT,dTdq[nq],dmudq[nq],
    dumudq[nq],dvmudq[nq],dkdq[nq];
  for (int n=0; n<npts; n++){
    iq       = n*nq;
    iqa      = n*nqa;
    iA       = n*nq*nq;
    r        = 1./q[iq];
    e        = q[iq+3]*r;
    p        = qa[iqa  ];
    u        = qa[iqa+1];
    v        = qa[iqa+2];
    t        = qa[iqa+3];
    mu       = qa[iqa+4];
    ui       = qai[iqa+1];
    vi       = qai[iqa+2];
    Ti       = qai[iqa+3];
    a        = xs[n]*xs[n]+ys[n]*ys[n];
    b        =(a+ys[n]*ys[n]/3.)/jac[n];
    c        =( -xs[n]*ys[n]/3.)/jac[n];
    d        =(a+xs[n]*xs[n]/3.)/jac[n];
    a        = a/jac[n];

    transport.getDViscosityDT   (1,&p,&t,&dmudT);
    transport.getDConductivityDT(1,&p,&t,&dkdT );

    f        = r*gm1/rGas;
    dTdq[0]  =-f*(e-(u*u+v*v));
    dTdq[1]  =-f*u;
    dTdq[2]  =-f*v;
    dTdq[3]  = f;

    for (int k=0; k<nq; k++) dmudq[k] = dmudT*dTdq[k];

    dumudq[0] = u*dmudq[0]-r*u*mu;
    dumudq[1] = u*dmudq[1]+r  *mu;
    dumudq[2] = u*dmudq[2]       ;
    dumudq[3] = u*dmudq[3]       ;

    dvmudq[0] = v*dmudq[0]-r*v*mu;
    dvmudq[1] = v*dmudq[1]       ;
    dvmudq[2] = v*dmudq[2]+r  *mu;
    dvmudq[3] = v*dmudq[3]       ;

    for (int k=0; k<nq; k++) dkdq[k] = dkdT*dTdq[k];

    A[iA   ] = 0.;
    A[iA+ 1] = 0.;
    A[iA+ 2] = 0.;
    A[iA+ 3] = 0.;

    A[iA+ 4] = dmudq[0]*(b*ui+c*vi);
    A[iA+ 5] = dmudq[1]*(b*ui+c*vi);
    A[iA+ 6] = dmudq[2]*(b*ui+c*vi);
    A[iA+ 7] = dmudq[3]*(b*ui+c*vi);

    A[iA+ 8] = dmudq[0]*(c*ui+d*vi);
    A[iA+ 9] = dmudq[1]*(c*ui+d*vi);
    A[iA+10] = dmudq[2]*(c*ui+d*vi);
    A[iA+11] = dmudq[3]*(c*ui+d*vi);

    A[iA+12] = dumudq[0]*(b*ui+c*vi)+dvmudq[0]*(c*ui+d*vi)+dkdq[0]*a*Ti;
    A[iA+13] = dumudq[1]*(b*ui+c*vi)+dvmudq[1]*(c*ui+d*vi)+dkdq[1]*a*Ti;
    A[iA+14] = dumudq[2]*(b*ui+c*vi)+dvmudq[2]*(c*ui+d*vi)+dkdq[2]*a*Ti;
    A[iA+15] = dumudq[3]*(b*ui+c*vi)+dvmudq[3]*(c*ui+d*vi)+dkdq[3]*a*Ti;
  }
}
