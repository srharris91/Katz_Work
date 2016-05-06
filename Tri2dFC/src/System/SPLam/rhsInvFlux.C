#include "Tri2dFCSPLam.h"


void Tri2dFCSPLam::rhsInvFlux(const int& npts,
			      const double& aL,
			      const double& aR,
			      const double* Ax,
			      const double* Ay,
			      const double* qL,
			      const double* qR,
			      const double* qaL,
			      const double* qaR,
			      double* f)
{
  int iq,iqa;
  double qn;
  for (int n=0; n<npts; n++){
    iq       = nq *n;
    iqa      = nqa*n;

    qn       = Ax[n]*qaL[iqa+1]+Ay[n]*qaL[iqa+2];
    f[iq  ]  = aL*(qn*qL[iq  ]               );
    f[iq+1]  = aL*(qn*qL[iq+1]+Ax[n]*qaL[iqa]);
    f[iq+2]  = aL*(qn*qL[iq+2]+Ay[n]*qaL[iqa]);
    f[iq+3]  = aL*(qn*qL[iq+3]+qn   *qaL[iqa]);

    qn       = Ax[n]*qaR[iqa+1]+Ay[n]*qaR[iqa+2];
    f[iq  ] += aR*(qn*qR[iq  ]               );
    f[iq+1] += aR*(qn*qR[iq+1]+Ax[n]*qaR[iqa]);
    f[iq+2] += aR*(qn*qR[iq+2]+Ay[n]*qaR[iqa]);
    f[iq+3] += aR*(qn*qR[iq+3]+qn   *qaR[iqa]);
  }
}
