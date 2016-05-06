#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::rhsTime(const int& j)
{
  int ne,nm,ni;
  double ns,jac1,s1;
  for (int n=0; n<nSurfNode; n++)
    for (int i=psp2S(n); i<psp2S(n+1); i++){
      ne   = psp1S(i,0);
      nm   = psp1S(i,1);
      ni   = surfElem(ne,nm);
      ns   = wsp1S(i);
      jac1 = jac(ne,nm,j);
      for (int k=0; k<nq; k++){
	s1        = bdf(0)*q(ni,j,k)*jac1;
	for (int m=0; m<timeAcc; m++) s1 += bdf(m+1)*qt(ne,nm,j,k,m);
	s1       /= dtUnsteady;
	r(n,j,k) += ns*s1;
      }}
}


void Strand2dFCBlockSolver::rhsTime()
{
  int ne,nm,ni;
  double ns,jac1,s1;
  for (int n=0; n<nSurfNode; n++)
    for (int i=psp2S(n); i<psp2S(n+1); i++){
      ne = psp1S(i,0);
      nm = psp1S(i,1);
      ni = surfElem(ne,nm);
      ns = wsp1S(i);
      for (int j=0; j<nStrandNode; j++){
	jac1 = jac(ne,nm,j);
	for (int k=0; k<nq; k++){
	  s1        = bdf(0)*q(ni,j,k)*jac1;
	  for (int m=0; m<timeAcc; m++) s1 += bdf(m+1)*qt(ne,nm,j,k,m);
	  s1       /= dtUnsteady;
	  r(n,j,k) += ns*s1;
	}}}
}
