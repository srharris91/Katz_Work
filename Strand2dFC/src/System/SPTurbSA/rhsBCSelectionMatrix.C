#include "Strand2dFCSPTurbSA.h"


void Strand2dFCSPTurbSA::rhsBCSelectionMatrix(const int& npts,
					const int* tag,
					const double* nx,
					const double* q,
					const double* qa,
					double* L)
{
  int ix,iq,iqa,iL;
  for (int n=0; n<npts; n++){
    ix  = ndim*n;
    iq  = nq*n;
    iqa = nqa*n;
    iL  = nq*nq*n;
    bc[tag[n]]->BCSelectionMatrix(nx+ix,
				  q +iq,
				  qa+iqa,
				  L +iL);
  }
}
