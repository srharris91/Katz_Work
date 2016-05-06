#include "Strand2dFCSPLam.h"


void Strand2dFCSPLam::lhsBCVectorSelfJacobian(const int& npts,
					   const int* tag,
					   const double* nx,
					   const double* q,
					   const double* qa,
					   double* M)
{
  int ix,iq,iqa,iM;
  for (int n=0; n<npts; n++){
    ix  = ndim*n;
    iq  = nq*n;
    iqa = nqa*n;
    iM  = nq*nq*n;
    bc[tag[n]]->BCVectorSelfJacobian(nx+ix,
				     q +iq,
				     qa+iqa,
				     M +iM);
  }
}
