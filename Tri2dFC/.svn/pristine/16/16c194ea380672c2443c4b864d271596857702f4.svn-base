#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::specRadv()
{
  // initialize spectral radius to zero
  radv.set(0.);


  // interior edge contributions
  int n1,n2;
  double Ax,Ay,sr,a=.5,b=.5;
  for(int n=0; n<nEdge; n++){
    n1        = edge(n,0);
    n2        = edge(n,1);
    Ax        = area(n,0);
    Ay        = area(n,1);
    sys->stepVisEigenvalue(1,a,b,&Ax,&Ay,&q(n1,0),&q(n2,0),
			   &qa(n1,0),&qa(n2,0),&v(n1),&v(n2),&sr);
    radv(n1) += sr;
    radv(n2) += sr;
  }

  // boundary edge contributions
  a = 5./6.;
  b = 1./6.;
  int m=0;
  for(int n=nEdge-nEdgeBd; n<nEdge; n++){
    n1        = edge(n,0);
    n2        = edge(n,1);
    Ax        = areaBd(m  ,0);
    Ay        = areaBd(m++,1);
    sys->stepVisEigenvalue(1,a,b,&Ax,&Ay,&q(n1,0),&q(n2,0),
			   &qa(n1,0),&qa(n2,0),&v(n1),&v(n2),&sr);
    radv(n1) += sr;
    sys->stepVisEigenvalue(1,b,a,&Ax,&Ay,&q(n1,0),&q(n2,0),
			   &qa(n1,0),&qa(n2,0),&v(n1),&v(n2),&sr);
    radv(n2) += sr;
  }
}
