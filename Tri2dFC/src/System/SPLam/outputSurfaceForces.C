#include "Tri2dFCSPLam.h"


void Tri2dFCSPLam::outputSurfaceForces(const int& npts,
				       const int* tag,
				       const double* A,
				       const double* q,
				       const double* qx,
				       const double* qy,
				       const double* qa,
				       const double* qax,
				       const double* qay,
				       double* force)
{
  int iA,iq,iqa,ifc;
  for (int n=0; n<npts; n++){
    iA   = ndim*n;
    iq   = nq*n;
    iqa  = nqa*n;
    ifc  = ndim*n;
    bc[tag[n]]->surfaceForces(A+iA,
			      q+iq,
			      qx+iq,
			      qy+iq,
			      qa+iqa,
			      qax+iqa,
			      qay+iqa,
			      force+ifc);
  }
}
