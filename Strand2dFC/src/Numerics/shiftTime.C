#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::shiftTime(const int& step)
{
  //NOTE: this assumes constant Jacobians in time, which will need to
  //be fixed for deforming grids
  int ni;
  double jaci;
  if (step == 1)
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<meshOrder+1; i++){
	ni = surfElem(n,i);
	for (int j=0; j<nStrandNode; j++){
	  jaci = jac(n,i,j);
	  for (int k=0; k<nq; k++)
	    for (int t=0; t<timeAcc; t++)
	      qt(n,i,j,k,t) = q(ni,j,k)*jaci;
	}}
  else{
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<meshOrder+1; i++)
	for (int j=0; j<nStrandNode; j++)
	  for (int k=0; k<nq; k++)
	    for (int t=timeAcc-1; t>0; t--)
	      qt(n,i,j,k,t) = qt(n,i,j,k,t-1);
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<meshOrder+1; i++){
	ni = surfElem(n,i);
	for (int j=0; j<nStrandNode; j++){
	  jaci = jac(n,i,j);
	  for (int k=0; k<nq; k++)
	    qt(n,i,j,k,0) = q(ni,j,k)*jaci;
	}}}
}
