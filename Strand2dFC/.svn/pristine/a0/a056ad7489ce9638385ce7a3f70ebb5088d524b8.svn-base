#include "Strand2dFCSPTurbSA.h"


void Strand2dFCSPTurbSA::initPenaltyData(const int& npts,
				         const int* tag,
				         const double* x,
				         double* f)
{
  int ix,iq;
  for (int n=0; n<npts; n++){
    ix  = ndim*n;
    iq  = nq*n;
    bc[tag[n]]->penaltyData(x+ix,
			    f+iq);
  }
}
