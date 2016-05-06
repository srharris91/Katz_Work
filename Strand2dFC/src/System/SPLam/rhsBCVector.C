#include "Strand2dFCSPLam.h"


void Strand2dFCSPLam::rhsBCVector(const int& npts,
				  const int* tag,
				  const double* nx,
				  const double* q,
				  const double* qa,
				  double* rb)
{
  int ix,iq,iqa;
  for (int n=0; n<npts; n++){
    ix  = ndim*n;
    iq  = nq*n;
    iqa = nqa*n;
    bc[tag[n]]->BCVector(nx+ix,
                         q +iq,
                         qa+iqa,
                         rb+iq);
  }
}
