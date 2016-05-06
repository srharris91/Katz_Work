#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::rhsBoundary(const int& j)
{
  if (j == 0) rhsBoundary();
}


void Strand2dFCBlockSolver::rhsBoundary()
{
  // compute update for strand root nodes (at this time, I treat viscous
  // wall strongly, which requires this older treatment)
  int j=0; // strand roots
  double rb[nq],L[nq*nq],a[nq];
  for (int n=0; n<nSurfNode; n++){
    sys->rhsBCVector(1,&surfNodeTag(n,0),&sn(n,0,0),
		     &q(n,j,0),&qa(n,j,0),&rb[0]);
    sys->rhsBCSelectionMatrix(1,&surfNodeTag(n,0),&sn(n,0,0),
			      &q(n,j,0),&qa(n,j,0),&L[0]);
    matmul(nq,nq,1,&L[0],&r(n,j,0),&a[0]);
    for (int k=0; k<nq; k++) r(n,j,k) = rb[k]+a[k];
  }
}
