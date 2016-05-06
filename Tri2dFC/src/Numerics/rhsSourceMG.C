#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::rhsSourceMG(const int& mode,
				     const int& stage)
{
  // compute forcing term
  if (stage == 0 && (mode == 2 || mode == 3)){
    double a,b;
    for(int n=0; n<nNode; n++)
      for (int k=0; k<nq; k++){
	a = 1.;
	if (fwc(n,k) < 0.) a = -1.;
	b = min(fabs(fwc(n,k)),fabs(r(n,k)));
	fwc(n,k) = relax*a*b-r(n,k);
	//fwc(n,k) = relax*fwc(n,k)-r(n,k);
      }
  }


  // add forcing term
  for(int n=0; n<nNode; n++)
    for (int k=0; k<nq; k++) r(n,k) += fwc(n,k);
}
