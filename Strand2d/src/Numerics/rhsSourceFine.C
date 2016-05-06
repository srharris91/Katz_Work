#include "StrandBlockSolver.h"


void StrandBlockSolver::rhsSourceFine()
{
  int npts=1;
  double f[nq];


  // compute source terms
  for (int n=0; n<nFaces-nGfaces; n++)
  for (int j=1; j<fClip(n)+1; j++){
    sys->rhsSource(npts,&q(0,j,n),&qa(0,j,n),&qx(0,0,j,n),&qax(0,0,j,n),&f[0]);
    for (int k=0; k<nq; k++) r(k,j,n) -= f[k]*v(j,n);
  }
}
