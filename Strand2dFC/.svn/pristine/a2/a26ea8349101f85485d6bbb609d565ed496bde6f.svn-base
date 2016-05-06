#include "Strand2dFCSPTurbSA.h"


void Strand2dFCSPTurbSA::stepVisEigenvalue(const int& npts,
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
  double mu,nu,k,AA,rho,Cpn,srL,srR,chi3,chi,fv1,rnu,mut;  
  
  for (int n=0; n<npts; n++){
    iq    = nq  *n;
    iqa   = nqa *n;
    Cpn   = gamma*rGas/gm1;
      
    rho  = qL [iq  ];
    rnu  = qL [iq+4];
    nu   = qaL[iqa+4];
    mu   = qaL[iqa+5];
    k    = qaL[iqa+6];
    chi  = rnu/mu;
    chi3 = chi*chi*chi;
    fv1  = chi3/(chi3+cv1*cv1*cv1);
    mut  = rnu*fv1;
    if (rnu < 0.) mut = 0.;
    AA    =(Ax[n]*Ax[n]+Ay[n]*Ay[n])/vL[n];
    srL   = AA*(mu+gamma*k/Cpn+mut*(1.+gamma)/PrnT)/qL[iq];

    rho  = qR [iq  ];
    rnu  = qR [iq+4];
    nu   = qaR[iqa+4];
    mu   = qaR[iqa+5];
    k    = qaR[iqa+6];
    chi  = rnu/mu;
    chi3 = chi*chi*chi;
    fv1  = chi3/(chi3+cv1*cv1*cv1);
    mut  = rnu*fv1;
    if (rnu < 0.) mut = 0.;
    AA    =(Ax[n]*Ax[n]+Ay[n]*Ay[n])/vR[n];
    srR   = AA*(mu+gamma*k/Cpn+mut*(1.+gamma)/PrnT)/qR[iq];

    sr[n] = aL*srL+aR*srR;
  }
}
