#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::rhsSourceMG(const int& mode,
					const int& stage,
					const int& sweep,
					const int& j)
{
  if (stage == 0 && sweep == 0 && (mode == 2 || mode == 3)){
    double a,b;
    for(int n=0; n<nSurfNode; n++)
      for (int k=0; k<nq; k++){
	a = 1.;
	if (fwc(n,j,k) < 0.) a = -1.;
	b = min(fabs(fwc(n,j,k)),fabs(r(n,j,k)));
	fwc(n,j,k) = relax*a*b-r(n,j,k);
      }}


  for(int n=0; n<nSurfNode; n++)
    for (int k=0; k<nq; k++) r(n,j,k) += fwc(n,j,k);
}


void Strand2dFCBlockSolver::rhsSourceMG()
{
  for(int n=0; n<nSurfNode; n++)
    for(int j=0; j<nStrandNode; j++)
      for (int k=0; k<nq; k++) r(n,j,k) += fwc(n,j,k);
}
