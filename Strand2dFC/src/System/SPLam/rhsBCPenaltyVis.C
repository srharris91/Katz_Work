#include "Strand2dFCSPLam.h"


void Strand2dFCSPLam::rhsBCPenaltyVis(const int& npts,
				      const int* tag,
				      const double* Pinv0,
				      const double* q,
				      const double* qa,
				      const double* gv,
				      const double* uw,
				      double* rb)
{
  int ix,iq,iqa;
  for (int n=0; n<npts; n++){
    ix  = ndim*n;
    iq  = nq*n;
    iqa = nqa*n;
    bc[tag[n]]->BCPenaltyVis(Pinv0[n],
			     q +iq,
			     qa+iqa,
			     gv+iq,
			     uw+ix,
			     rb+iq);
  }
}
