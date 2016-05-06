#include "Strand2dFCSPLam.h"


void Strand2dFCSPLam::outputSurfaceForces(const int& npts,
					  const int* tag,
					  const double* xs,
					  const double* ys,
					  const double* q,
					  const double* qa,
					  const double* qx,
					  const double* qy,
					  const double* qax,
					  const double* qay,
					  double* force)
{
  int iq,iqa,ifc;
  for (int n=0; n<npts; n++){
    iq   = nq*n;
    iqa  = nqa*n;
    ifc  = ndim*n;
    bc[tag[n]]->surfaceForces(xs+n,
			      ys+n,
			      q+iq,
			      qa+iqa,
			      qx+iq,
			      qy+iq,
			      qax+iqa,
			      qay+iqa,
			      force+ifc);
  }
}
