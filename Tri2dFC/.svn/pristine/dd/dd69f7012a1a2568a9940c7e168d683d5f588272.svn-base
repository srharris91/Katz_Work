#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::shiftTime(const int& step)
{
  if (step == 1)
    for (int t=0; t<timeAcc; t++)
      for (int n=0; n<nNode; n++)
	for (int k=0; k<nq; k++)
	  qt(t,n,k) = q(n,k);
  else{
    for (int t=timeAcc-1; t>0; t--)
      for (int n=0; n<nNode; n++)
	for (int k=0; k<nq; k++)
	  qt(t,n,k) = qt(t-1,n,k);
    for (int n=0; n<nNode; n++)
      for (int k=0; k<nq; k++)
	qt(0,n,k) = q(n,k);
  }
}
