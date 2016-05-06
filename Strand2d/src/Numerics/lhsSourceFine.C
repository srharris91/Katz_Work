#include "StrandBlockSolver.h"


void StrandBlockSolver::lhsSourceFine()
{
  // compute source terms
  int npts=1,m;
  double a[nq*nq];
  for (int n=0; n<nFaces-nGfaces; n++)
  for (int j=1; j<fClip(n)+1; j++){
    sys->lhsSourceJacobian(npts,&v(j,n),&q(0,j,n),&qa(0,j,n),&qx(0,0,j,n),
			   &qax(0,0,j,n),&a[0]);
    for (int k=0; k<nq; k++){
      m = k*nq;
      for (int l=0; l<nq; l++) dd(l,k,j,n) -=(a[m+l]*v(j,n));
    }}
}
