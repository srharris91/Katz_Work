#include "StrandSPTurbSA.h"


void StrandSPTurbSA::rhsBCVector(const int& npts,
				 const int* tag,
				 const double* nx,
				 const double* wx,
				 const double* qe,
				 const double* qae,
				 const double* q,
				 const double* qa,
				 double* r)
{
  int ix,iq,iqa;
  for (int n=0; n<npts; n++){
    ix   = ndim*n;
    iq   = nq*n;
    iqa  = nqa*n;
    bc[tag[n]]->BCVector(nx +ix,
			 wx +ix,
			 qe +iq,
			 qae+iqa,
			 q  +iq,
			 qa +iqa,
			 r  +iq);
  }
}
