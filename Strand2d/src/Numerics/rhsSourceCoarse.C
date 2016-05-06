#include "StrandBlockSolver.h"


void StrandBlockSolver::rhsSourceCoarse()
{
  int npts=1;
  double f[nq];


  // compute source terms
  for (int n=0; n<nFaces-nGfaces; n++)
  for (int j=1; j<fClip(n)+1; j++){
    sys->rhsSourceCoarse(npts,&v(j,n),&q(0,j,n),&qa(0,j,n),&f[0]);
    for (int k=0; k<nq; k++) r(k,j,n) -= f[k]*v(j,n);
  }
}
