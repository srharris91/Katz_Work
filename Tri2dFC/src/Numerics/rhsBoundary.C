#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::rhsBoundary()
{
  int m=0;
  double rb[nq],L[nq*nq],a[nq];
  for (int n=nNode-nNodeBd; n<nNode; n++){
    sys->rhsBCVector(1,&nodeBd(m),&ln(m,0),&q(n,0),&qa(n,0),&rb[0]);
    sys->rhsBCSelectionMatrix(1,&nodeBd(m),&ln(m,0),&q(n,0),&qa(n,0),&L[0]);
    matmul(nq,nq,1,&L[0],&r(n,0),&a[0]);
    for (int k=0; k<nq; k++) r(n,k) = rb[k]+a[k];
    m++;
  }
}
