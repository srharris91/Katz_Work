#include "StrandSPTurbSA.h"


void StrandSPTurbSA::rhsVisFluxCoarse(const int& npts,
				   const double* A,
				   const double* vl,
				   const double* vr,
				   const double* ql,
				   const double* qr,
				   const double* qal,
				   const double* qar,
				   double* f)
{
  int iq,iqa,iM,j=1;
  double qe[nq*npts],qae[nqa*npts],ve[npts],M[nq*nq*npts],dq[nq];
  for (int n=0; n<npts; n++){
    iq    = nq *n;
    iqa   = nqa*n;
    for (int k=0; k<nq;  k++) qe [iq +k] = .5*(ql [iq +k]+qr [iq +k]);
    for (int k=0; k<nqa; k++) qae[iqa+k] = .5*(qal[iqa+k]+qar[iqa+k]);
    ve[n] = .5*(vl[n]+vr[n]);
  }
   
  lhsVisFluxJacobian(npts,A,A,&qe[0],&qae[0],&M[0]);

  for (int n=0; n<npts; n++){
    iq    = nq *n;
    iqa   = nqa*n;
    iM    = nq*nq*n;
    for (int k=0; k<nq; k++) dq[k] = (qar[iqa+k]-qal[iqa+k])/ve[n];
    matmul(nq,nq,j,&M[iM],&dq[0],&f[iq]);
  }
}
