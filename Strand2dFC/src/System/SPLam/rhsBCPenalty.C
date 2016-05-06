#include "Strand2dFCSPLam.h"


void Strand2dFCSPLam::rhsBCPenalty(const int& npts,
				   const int* tag,
				   const int& inout,
				   const double* A,
				   const double* Pinv0,
				   const double* q,
				   const double* qa,
				   const double* g,
				   const double* uw,
				   double* rb)
{
  int ix,iq,iqa;
  for (int n=0; n<npts; n++){
    ix  = ndim*n;
    iq  = nq*n;
    iqa = nqa*n;
    bc[tag[n]]->BCPenalty(inout,
			  A +ix,
			  Pinv0[n],
			  q +iq,
			  qa+iqa,
			  g +iq,
			  uw+ix,
			  rb+iq);
  }
}
