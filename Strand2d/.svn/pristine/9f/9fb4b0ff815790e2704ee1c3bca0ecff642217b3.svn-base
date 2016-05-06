#include "StrandSPLam.h"


void StrandSPLam::rhsInvFlux(const int& npts,
			     const double* A,
			     const double* xv,
			     const double* ql,
			     const double* qr,
			     double* f)
{
  // inviscid flux terms
  int iq,iA;
  double Ax,Ay,qsl,qsr,rhl,rhr,pl,pr;
  for (int n=0; n<npts; n++){
    iq      = nq  *n;
    iA      = ndim*n;
    Ax      = A[iA  ];
    Ay      = A[iA+1];
    pl      = gm1*(ql[iq+3]-.5*(ql[iq+1]*ql[iq+1]+ql[iq+2]*ql[iq+2])/ql[iq]);
    pr      = gm1*(qr[iq+3]-.5*(qr[iq+1]*qr[iq+1]+qr[iq+2]*qr[iq+2])/qr[iq]);
    qsl     =(Ax*ql[iq+1]+Ay*ql[iq+2])/ql[iq]-xv[n];
    qsr     =(Ax*qr[iq+1]+Ay*qr[iq+2])/qr[iq]-xv[n];
    rhl     = ql[iq+3]+pl;
    rhr     = qr[iq+3]+pr;
    f[iq  ] = qsl*ql[iq  ]+qsr*qr[iq  ]           ;
    f[iq+1] = qsl*ql[iq+1]+qsr*qr[iq+1]+Ax*(pl+pr);
    f[iq+2] = qsl*ql[iq+2]+qsr*qr[iq+2]+Ay*(pl+pr);
    f[iq+3] = qsl*rhl     +qsr*rhr                ;
    for (int k=0; k<nq; k++) f[iq+k] *= .5;
  }
}
