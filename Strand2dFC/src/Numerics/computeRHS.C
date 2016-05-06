#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::computeRHS(const int& step,
				       const int& pseudoStep)
{
  r.set(0.); // initialize residuals
  d.set(0.);

  if (step > 0) rhsTime(); // physical-time terms

  if (dissipation != 0) // dissipation terms
    for(int j=0; j<nStrandNode; j++) rhsDissipation(j);

  if (inviscid != 0) // inviscid terms
    for(int j=0; j<nStrandNode; j++) rhsInviscid(j);

  if (viscous != 0) // viscous terms
    for(int j=0; j<nStrandNode; j++) rhsViscous(j);

  if (source == 1) // physical source terms
    for(int j=0; j<nStrandNode; j++) rhsSource(j);

  if (sourceMMS == 1) rhsSourceMMS(); // MMS source terms

  //Add convective and diffusive residual vectors
  for(int n=0; n<nSurfNode; n++)
    for(int j=0; j<nStrandNode; j++)
      for (int k=0; k<nq; k++) r(n,j,k) += d(n,j,k);

  rhsBoundary(); // strong boundaries

  if (level > 0) rhsSourceMG(); // multigrid forcing terms
}
