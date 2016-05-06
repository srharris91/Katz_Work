#include "Strand2dFCSPTurbSA.h"


void Strand2dFCSPTurbSA::lhsBCPenaltyJacobian(const int& npts,
					      const int* tag,
					      const int& inout,
					      const double* A,
					      const double* Pinv0,
					      const double* q,
					      const double* qa,
					      const double* g,
					      const double* uw,
					      double* M)
{
  int ix,iq,iqa,iM;
  for (int n=0; n<npts; n++){
    ix  = ndim*n;
    iq  = nq*n;
    iqa = nqa*n;
    iM  = n*nq*nq;
    bc[tag[n]]->BCPenaltyJacobian(inout,
				  A +ix,
				  Pinv0[n],
				  q +iq,
				  qa+iqa,
				  g +iq,
				  uw+ix,
				  M+iM);
  }
}
