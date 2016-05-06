#include "StrandSPTurbSA.h"


void StrandSPTurbSA::outputSurfaceForces(const int& npts,
				      const int* tag,
				      const double* A,
				      const double* q,
				      const double* qa,
				      const double* qx,
				      const double* qax,
				      double* force)
{
  int iA,iq,iqa,iqx,iqax,ifc;
  for (int n=0; n<npts; n++){
    iA   = ndim*n;
    iq   = nq*n;
    iqa  = nqa*n;
    iqx  = nq*ndim*n;
    iqax = nqa*ndim*n;
    ifc  = ndim*n;
    bc[tag[n]]->surfaceForces(A+iA,
			      q+iq,
			      qa+iqa,
			      qx+iqx,
			      qax+iqax,
			      force+ifc);
  }
}
