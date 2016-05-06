#include "Strand2dFCBlockSolver.h"
#include "matinv.h"


void Strand2dFCBlockSolver::computeLHS(const int& step,
				       const int& pseudoStep,
				       const int& sweep,
				       const int& subStep,
				       const int& stage,
				       const int& mode,
				       const int& j)
{
  // initialize diagonal Jacobian terms to zero
  Ads.set(0.);


  // physical-time terms
  //if (step > 0) lhsTime();


  // dissipation terms
  if (dissipation != 0) lhsDissipation(j);


  // inviscid terms
  if (inviscid != 0) lhsInviscid(j);


  // viscous terms
  //if (viscous != 0) lhsViscous(j);


  // physical source terms
  //if (source == 1) lhsSource();


  // pseudo-time terms
  /*
  pseudoTime(j);
  double ap,au,b=0.,c,eps=1.e-14,cfln=100000000000000.*cfl;
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<meshOrder+1; i++){
      ap = 1./(dt(n,j)/cfl*cfln);
      au = b/max(dtUnsteady,eps);
      c  = v(n,j)*(ap+au);
      for (int k=0; k<nq; k++) Ads(n,i,k,k) += c;
    }
  */
}
