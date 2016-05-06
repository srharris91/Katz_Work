#include "Tri2dFCSPLam.h"


void Tri2dFCSPLam::stepInvEigenvalue(const int& npts,
				     const double& aL,
				     const double& aR,
				     const double* Ax,
				     const double* Ay,
				     const double* qL,
				     const double* qR,
				     const double* qaL,
				     const double* qaR,
				     double* sr)
{
  int iqa;
  double t,c,qn;
  for (int n=0; n<npts; n++){
    iqa   = nqa*n;
    t     = .5*(qaL[iqa+3]+qaR[iqa+3]);
    c     = sqrt(gamma*rGas*t);
    c    *= sqrt(Ax[n]*Ax[n]+Ay[n]*Ay[n]);
    qn    = .5*Ax[n]*(qaL[iqa+1]+qaR[iqa+1])
          + .5*Ay[n]*(qaL[iqa+2]+qaR[iqa+2]);
    sr[n] = fabs(qn)+c;
  }
}
