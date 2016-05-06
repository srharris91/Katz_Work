#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::computeRHS(const int& step,
				    const int& pseudoStep,
				    const int& stage,
				    const int& mode)
{
  gradQComputed   = false;
  gradQaComputed  = false;
  limiterComputed = false;

  double a = rkb(stage);
  double b = 1.-a;


  // initialize diffusive and convective residuals
  if (stage == 0) dn.set(0.);
  else
    for(int n=0; n<nNode; n++)
      for (int k=0; k<nq; k++) dn(n,k) = d(n,k);
  r.set(0.);
  d.set(0.);


  // physical-time terms
  if (step > 0) rhsTime();


  // dissipation terms
  if (dissipation != 0 && a != 0.) rhsDissipation();


  // inviscid terms
  if (inviscid != 0) rhsInviscid();


  // viscous terms
  if (viscous != 0 && a != 0.) rhsViscous();


  // physical source terms
  if (source == 1) rhsSource();


  // MMS source terms
  if (sourceMMS == 1) rhsSourceMMS();


  //Add convective and diffusive residual vectors
  for(int n=0; n<nNode; n++)
    for (int k=0; k<nq; k++){
      d(n,k)  = a*d(n,k)+b*dn(n,k);
      r(n,k) += d(n,k);
    }


  // modify residuals at boundaries to incorporate boundary conditions
  rhsBoundary();


  // multigrid forcing terms on coarse levels
  //if (level > 0) rhsSourceMG(mode,
  //			     stage);
}
