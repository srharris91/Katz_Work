#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::specRadi()
{
  // initialize spectral radius to zero
  radi.set(0.);


  // interior edge contributions
  int n1,n2;
  double Ax,Ay,sr,a=.5,b=.5;
  for(int n=0; n<nEdge; n++){
    n1        = edge(n,0);
    n2        = edge(n,1);
    Ax        = area(n,0);
    Ay        = area(n,1);
    sys->stepInvEigenvalue(1,a,b,&Ax,&Ay,&q(n1,0),&q(n2,0),&qa(n1,0),&qa(n2,0),&sr);
    radi(n1) += sr;
    radi(n2) += sr;
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
    sys->stepInvEigenvalue(1,a,b,&Ax,&Ay,&q(n1,0),&q(n2,0),&qa(n1,0),&qa(n2,0),&sr);
    radi(n1) += sr;
    sys->stepInvEigenvalue(1,b,a,&Ax,&Ay,&q(n1,0),&q(n2,0),&qa(n1,0),&qa(n2,0),&sr);
    radi(n2) += sr;
  }
}
