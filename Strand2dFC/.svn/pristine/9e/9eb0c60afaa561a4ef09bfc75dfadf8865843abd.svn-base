#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::rhsSourceMMS(const int& j)
{
  int ne,ni,nm;
  double ns;
  for (int n=0; n<nSurfNode; n++)
    for (int i=psp2S(n); i<psp2S(n+1); i++){
      ne = psp1S(i,0);
      ni = psp1S(i,1);
      nm = surfElem(ne,ni);
      ns = wsp1S(i);
      for (int k=0; k<nq; k++) r(n,j,k) -= ns*jac(ne,ni,j)*s(nm,j,k);
    }
}


void Strand2dFCBlockSolver::rhsSourceMMS()
{
  int ne,ni,nm;
  double ns;
  for (int n=0; n<nSurfNode; n++)
    for (int i=psp2S(n); i<psp2S(n+1); i++){
      ne = psp1S(i,0);
      ni = psp1S(i,1);
      nm = surfElem(ne,ni);
      ns = wsp1S(i);
      for(int j=0; j<nStrandNode; j++)
	for (int k=0; k<nq; k++) r(n,j,k) -= ns*jac(ne,ni,j)*s(nm,j,k);
    }
}
