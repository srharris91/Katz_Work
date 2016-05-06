#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::specRadv(const int& j)
{
  // initialize spectral radius to zero
  radv.set(0.);


  // surface edge contributions
  int i1,i2,n1,n2;
  double Ax,Ay,sr,a=.5,b=.5;
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<nElemEdge; i++){
      i1 = elemEdge(i,0);
      i2 = elemEdge(i,1);
      n1 = surfElem(n,i1);
      n2 = surfElem(n,i2);
      Ax = .5*(yn(n1,j)+yn(n2,j));
      Ay =-.5*(xn(n1,j)+xn(n2,j));
      sys->stepVisEigenvalue(1,a,b,&Ax,&Ay,&q(n1,j,0),&q(n2,j,0),
			     &qa(n1,j,0),&qa(n2,j,0),&v(n1,j),&v(n2,j),&sr);
      radv(n1) += sr;
      radv(n2) += sr;
    }

  for (int n=0; n<nBndNode; n++){
    n1 = bndNode(n);
    Ax = yn(n1,j);
    Ay =-xn(n1,j);
    sys->stepVisEigenvalue(1,a,b,&Ax,&Ay,&q(n1,j,0),&q(n1,j,0),
			   &qa(n1,j,0),&qa(n1,j,0),&v(n1,j),&v(n1,j),&sr);
    radv(n1) += sr;
  }
}
