#include "Tri2dFCSPLam.h"


void Tri2dFCSPLam::outputSurfaceSolution(const int& npts,
					 ofstream& ffile,
					 const int* tag,
					 const double* x,
					 const double* y,
					 const double* q,
					 const double* qa,
					 const double* qx,
					 const double* qy,
					 const double* qax,
					 const double* qay)
{
  int ix,iy,iq,iqa;
  for (int n=0; n<npts; n++){
    ix   = n;
    iy   = n;
    iq   = nq*n;
    iqa  = nqa*n;
    bc[tag[n]]->surfaceSolution(ffile,
				x+ix,
				y+iy,
				q+iq,
				qa+iqa,
				qx+iq,
				qy+iq,
				qax+iqa,
				qay+iqa);
  }
}
