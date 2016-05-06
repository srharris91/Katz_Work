#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::specRadi(const int& j)
{
  // initialize spectral radius to zero
  radi.set(0.);


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
      sys->stepInvEigenvalue(1,a,b,&Ax,&Ay,&q(n1,j,0),&q(n2,j,0),
			     &qa(n1,j,0),&qa(n2,j,0),&sr);
      radi(n1) += sr;
      radi(n2) += sr;
    }

  for (int n=0; n<nBndNode; n++){
    n1 = bndNode(n);
    Ax = yn(n1,j);
    Ay =-xn(n1,j);
    sys->stepInvEigenvalue(1,a,b,&Ax,&Ay,&q(n1,j,0),&q(n1,j,0),
			   &qa(n1,j,0),&qa(n1,j,0),&sr);
    radi(n1) += sr;
  }


  /*
  // strand edge contributions
  if (j == 0){
    double dnr=1./deltaN;
    for (int n=0; n<nSurfNode; n++){
      Ax =-dnr*ysA(n,j);
      Ay = dnr*xsA(n,j);
      sys->stepInvEigenvalue(1,a,b,&Ax,&Ay,&q(n,j,0),&q(n,j,0),
			     &qa(n,j,0),&qa(n,j,0),&sr);
      radi(n) += sr;
      
      Ax =-dnr*(ysA(n,j)+ysA(n,j+1));
      Ay = dnr*(xsA(n,j)+xsA(n,j+1));
      sys->stepInvEigenvalue(1,a,b,&Ax,&Ay,&q(n,j,0),&q(n,j+1,0),
			     &qa(n,j,0),&qa(n,j+1,0),&sr);
      radi(n) += sr;
    }}

  else if (j == nStrandNode-1){
    double dnr=1./deltaN;
    for (int n=0; n<nSurfNode; n++){
      Ax =-dnr*(ysA(n,j-1)+ysA(n,j));
      Ay = dnr*(xsA(n,j-1)+xsA(n,j));
      sys->stepInvEigenvalue(1,a,b,&Ax,&Ay,&q(n,j-1,0),&q(n,j,0),
			     &qa(n,j-1,0),&qa(n,j,0),&sr);
      radi(n) += sr;
      
      Ax =-dnr*ysA(n,j);
      Ay = dnr*xsA(n,j);
      sys->stepInvEigenvalue(1,a,b,&Ax,&Ay,&q(n,j,0),&q(n,j,0),
			     &qa(n,j,0),&qa(n,j,0),&sr);
      radi(n) += sr;
    }}

  else{
    double dnr=.5/deltaN;
    for (int n=0; n<nSurfNode; n++){
      Ax =-dnr*(ysA(n,j-1)+ysA(n,j));
      Ay = dnr*(xsA(n,j-1)+xsA(n,j));
      sys->stepInvEigenvalue(1,a,b,&Ax,&Ay,&q(n,j-1,0),&q(n,j,0),
			     &qa(n,j-1,0),&qa(n,j,0),&sr);
      radi(n) += sr;
      
      Ax =-dnr*(ysA(n,j)+ysA(n,j+1));
      Ay = dnr*(xsA(n,j)+xsA(n,j+1));
      sys->stepInvEigenvalue(1,a,b,&Ax,&Ay,&q(n,j,0),&q(n,j+1,0),
			     &qa(n,j,0),&qa(n,j+1,0),&sr);
      radi(n) += sr;
    }}
  */
}
