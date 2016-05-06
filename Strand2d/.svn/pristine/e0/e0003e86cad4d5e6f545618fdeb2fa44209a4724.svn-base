#include "StrandBlockSolver.h"


void StrandBlockSolver::rhsTimeFine()
{
  double t0,t1,t2,t0v,t1v,t2v;
  t0 = 1.5/dtUnsteady;
  t1 =-2.0/dtUnsteady;
  t2 = 0.5/dtUnsteady;
  for (int n=0; n<nFaces-nGfaces; n++){
  for (int j=1; j<fClip(n)+1; j++){
    t0v       = t0*v (j,n);
    t1v       = t1*v1(j,n);
    t2v       = t2*v2(j,n);
    for (int k=0; k<nq; k++)
      r(k,j,n) += (t0v*q(k,j,n)+t1v*q1(k,j,n)+t2v*q2(k,j,n));
  }}
}
