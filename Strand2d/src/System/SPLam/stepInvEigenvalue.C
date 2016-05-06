#include "StrandSPLam.h"


void StrandSPLam::stepInvEigenvalue(const int& npts,
				    const double* A,
				    const double* xv,
				    const double* q1,
				    const double* q2,
				    const double* qa1,
				    const double* qa2,
				    double* sr)
{
  int iqa,iA,j=1;
  double p,t,c,Ax,Ay,qs;
  for (int n=0; n<npts; n++){
    iqa   = nqa *n;
    iA    = ndim*n;
    p     = .5*(qa1[iqa  ]+qa2[iqa  ]);
    t     = .5*(qa1[iqa+3]+qa2[iqa+3]);
    c     = sqrt(gamma*rGas*t);
    Ax    = A[iA  ];
    Ay    = A[iA+1];
    c    *= sqrt(Ax*Ax+Ay*Ay);
    qs    = .5*(Ax*(qa1[iqa+1]+qa2[iqa+1])+
	      Ay*(qa1[iqa+2]+qa2[iqa+2]));
    sr[n] = fabs(qs-xv[n])+c;
  }
}
