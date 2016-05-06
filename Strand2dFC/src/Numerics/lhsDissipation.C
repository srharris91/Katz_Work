#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::lhsDissipation(const int& j)
{
  // strand edge dissipation contributions
  int ni;
  double Ax,Ay,dnr=.5/deltaN;
  Array2D<double> A(nq,nq);

  if (j == 0) //root edge
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<meshOrder+1; i++){
	ni = surfElem(n,i);
	Ax =-.5*(ys(n,i,j)+ys(n,i,j+1));
	Ay = .5*(xs(n,i,j)+xs(n,i,j+1));
	sys->lhsDisFluxJacobian(1,&Ax,&Ay,&q(ni,j,0),&q(ni,j+1,0),&A(0,0));
	for (int k=0; k<nq; k++)
	  for (int l=0; l<nq; l++) Ads(n,i,k,l) += 2.*dnr*A(k,l);
      }

  else if (j == nStrandNode-1) //tip edge
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<meshOrder+1; i++){
	ni = surfElem(n,i);
	Ax =-.5*(ys(n,i,j-1)+ys(n,i,j));
	Ay = .5*(xs(n,i,j-1)+xs(n,i,j));
	sys->lhsDisFluxJacobian(1,&Ax,&Ay,&q(ni,j-1,0),&q(ni,j,0),&A(0,0));
	for (int k=0; k<nq; k++)
	  for (int l=0; l<nq; l++) Ads(n,i,k,l) += 2.*dnr*A(k,l);
      }

  else //interior edges
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<meshOrder+1; i++){
	ni = surfElem(n,i);
	Ax =-.5*(ys(n,i,j-1)+ys(n,i,j));
	Ay = .5*(xs(n,i,j-1)+xs(n,i,j));
	sys->lhsDisFluxJacobian(1,&Ax,&Ay,&q(ni,j-1,0),&q(ni,j,0),&A(0,0));
	for (int k=0; k<nq; k++)
	  for (int l=0; l<nq; l++) Ads(n,i,k,l) += dnr*A(k,l);
	Ax =-.5*(ys(n,i,j)+ys(n,i,j+1));
	Ay = .5*(xs(n,i,j)+xs(n,i,j+1));
	sys->lhsDisFluxJacobian(1,&Ax,&Ay,&q(ni,j,0),&q(ni,j+1,0),&A(0,0));
	for (int k=0; k<nq; k++)
	  for (int l=0; l<nq; l++) Ads(n,i,k,l) += dnr*A(k,l);
      }


  // clean up
  A.deallocate();
}
