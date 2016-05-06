#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::lhsInviscid(const int& j)
{
  int ni;
  double uw[2],A[2],pb[nq],dnr=1./deltaN,Pdnr=Pinv0*dnr;
  Array2D<double> Ab(nq,nq);
  uw[0] = 0.; //for now, use fixed walls in boundary conditions
  uw[1] = 0.;

  // strand boundary edge inviscid flux contributions
  if (j == 0) // strand root
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<meshOrder+1; i++){
	ni = surfElem(n,i);

	// boundary flux Jacobian
	A[0] =-ys(n,i,j);
	A[1] = xs(n,i,j);
	sys->lhsInvFluxJacobian(1,&A[0],&A[1],&q(ni,j,0),&qa(ni,j,0),&Ab(0,0));
	for (int k=0; k<nq; k++)
	  for (int l=0; l<nq; l++) Ads(n,i,k,l) -= dnr*Ab(k,l);

	// boundary penalty
	sys->lhsBCPenaltyJacobian(1,&surfNodeTag(ni,0),1,&A[0],&Pdnr,
				  &q(ni,j,0),&qa(ni,j,0),&surfData(ni,0,0),
				  &uw[0],&Ab(0,0));
	for (int k=0; k<nq; k++)
	  for (int l=0; l<nq; l++) Ads(n,i,k,l) -= Pdnr*Ab(k,l);
      }

  if (j == nStrandNode-1) // strand tip
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<meshOrder+1; i++){
	ni = surfElem(n,i);

	// boundary flux Jacobian
	A[0] =-ys(n,i,j);
	A[1] = xs(n,i,j);
	sys->lhsInvFluxJacobian(1,&A[0],&A[1],&q(ni,j,0),&qa(ni,j,0),&Ab(0,0));
	for (int k=0; k<nq; k++)
	  for (int l=0; l<nq; l++) Ads(n,i,k,l) += dnr*Ab(k,l);

	// boundary penalty
	sys->lhsBCPenaltyJacobian(1,&surfNodeTag(ni,1),-1,&A[0],&Pdnr,
				  &q(ni,j,0),&qa(ni,j,0),&surfData(ni,1,0),
				  &uw[0],&Ab(0,0));
	for (int k=0; k<nq; k++)
	  for (int l=0; l<nq; l++) Ads(n,i,k,l) += Pdnr*Ab(k,l);
      }


  // clean up
  Ab.deallocate();
}
