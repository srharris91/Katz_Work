#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::solve(const int& step,
				  const int& pseudoStep,
				  const int& sweep,
				  const int& stage,
				  const int& mode,
				  const int& j)
{
  // compute RHS terms
  double a=rkb(stage),b=1.-a;

  // initialize diffusive and convective residuals
  if (stage == 0) for(int n=0; n<nSurfNode; n++)
		    for (int k=0; k<nq; k++) dn(n,j,k) = 0.;
  else for(int n=0; n<nSurfNode; n++)
	 for (int k=0; k<nq; k++) dn(n,j,k) = d(n,j,k);
  for(int n=0; n<nSurfNode; n++)
    for (int k=0; k<nq; k++) r(n,j,k) = 0.;
  for(int n=0; n<nSurfNode; n++)
    for (int k=0; k<nq; k++) d(n,j,k) = 0.;

  if (step > 0) rhsTime(j); // physical-time terms

  if (dissipation != 0 && a != 0.) rhsDissipation(j); // dissipation terms

  if (inviscid != 0) rhsInviscid(j); // inviscid terms

  if (viscous != 0 && a != 0.) rhsViscous(j); // viscous terms

  if (source == 1) rhsSource(j); // physical source terms

  if (sourceMMS == 1) rhsSourceMMS(j); // MMS source terms

  //Add convective and diffusive residual vectors
  for(int n=0; n<nSurfNode; n++)
    for (int k=0; k<nq; k++){
      d(n,j,k)  = a*d(n,j,k)+b*dn(n,j,k);
      r(n,j,k) += d(n,j,k);
    }

  rhsBoundary(j); // strong boundary conditions

  if (level > 0) rhsSourceMG(mode,stage,sweep,j); // multigrid forcing terms


  // compute LHS on first stage
  if (stage == 0){
    Ads.set(0.); // initialize diagonal Jacobian terms to zero

    if (dissipation != 0) lhsDissipation(j); // dissipation terms

    if (inviscid != 0) lhsInviscid(j); // inviscid terms

    if (viscous != 0) lhsViscous(j); // viscous terms

    if (source == 1) lhsSource(j); // physical source terms

    //lumped treatment of Jacobian matrix for LHS
    int i1,i2,n1,n2;
    double ds5=.5*deltaS;
    Adl.set(0.);
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<nElemEdge; i++){
	i1 = elemEdge(i,0);
	i2 = elemEdge(i,1);
	n1 = surfElem(n,i1);
	n2 = surfElem(n,i2);
	for (int k=0; k<nq; k++)
	  for (int l=0; l<nq; l++){
	    Adl(n1,k,l) += ds5*Ads(n,i1,k,l);
	    Adl(n2,k,l) += ds5*Ads(n,i2,k,l);
	  }}}


  // save first stage q, on the fine level save q for RMS purposes
  if (stage == 0){
    for(int n=0; n<nSurfNode; n++)
      for (int k=0; k<nq; k++) qn(n,j,k) = q(n,j,k);
    if (level == 0 && sweep == 0)
      for(int n=0; n<nSurfNode; n++)
	for (int k=0; k<nq; k++) q0(n,j,k) = q(n,j,k);
  }


  // update
  update(step,stage,j);
}
