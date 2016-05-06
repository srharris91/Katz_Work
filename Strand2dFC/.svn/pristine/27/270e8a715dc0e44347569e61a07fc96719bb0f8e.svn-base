#include "Strand2dFCSPLam.h"


void Strand2dFCSPLam::stepVisEigenvalue(const int& npts,
				     const double& aL,
				     const double& aR,
				     const double* Ax,
				     const double* Ay,
				     const double* qL,
				     const double* qR,
				     const double* qaL,
				     const double* qaR,
				     const double* vL,
				     const double* vR,
				     double* sr)
{
  int iq,iqa;
  double mu,k,AA,Cpn,srL,srR;
  for (int n=0; n<npts; n++){
    iq    = nq  *n;
    iqa   = nqa *n;
    Cpn   = gamma*rGas/gm1;

    mu    = qaL[iqa+4];
    k     = qaL[iqa+5];
    AA    =(Ax[n]*Ax[n]+Ay[n]*Ay[n])/vL[n];
    srL   = AA*(mu+gamma*k/Cpn)/qL[iq];

    mu    = qaR[iqa+4];
    k     = qaR[iqa+5];
    AA    =(Ax[n]*Ax[n]+Ay[n]*Ay[n])/vR[n];
    srR   = AA*(mu+gamma*k/Cpn)/qR[iq];

    sr[n] = aL*srL+aR*srR;
  }
}
