#include "StrandBlockSolver.h"


void StrandBlockSolver::rhsTimeCoarse()
{
  double t0,t0v;
  t0 = 1.5/dtUnsteady;
  for (int n=0; n<nFaces-nGfaces; n++){
  for (int j=1; j<fClip(n)+1; j++){
    t0v = t0*v(j,n);
    for (int k=0; k<nq; k++)
      r(k,j,n) += (t0v*q(k,j,n));
  }}
}
