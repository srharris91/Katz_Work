#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::rhsSourceMMS()
{
  if      (order == 1)
    for (int n=0; n<nNode; n++)
      for (int k=0; k<nq; k++) r(n,k) -= v(n)*s(n,k);

  else if (order == 2){
    int n1,n2;
    double Ax,Ay,dx,dy,vv,a;
    for (int n=0; n<nEdge; n++){
      n1 = edge(n,0);
      n2 = edge(n,1);
      Ax = area(n,0);
      Ay = area(n,1);
      dx = x(n2,0)-x(n1,0);
      dy = x(n2,1)-x(n1,1);
      vv = .125*(dx*Ax+dy*Ay);
      for (int k=0; k<nq; k++){
	a        = vv*(s(n1,k)+s(n2,k));
	r(n1,k) -= a;
	r(n2,k) -= a;
      }}}

  else if (order == 3){
    int n1,n2;
    double Ax,Ay,dx,dy,vv,ds0,ds1,ds2;
    Array3D<double> sxx(nNode,5,nq);

    // compute gradient and Hessian of S
    hessian(nq,&s(0,0),&sxx(0,0,0));

    for (int n=0; n<nEdge; n++){
      n1 = edge(n,0);
      n2 = edge(n,1);
      Ax = area(n,0);
      Ay = area(n,1);
      dx = x(n2,0)-x(n1,0);
      dy = x(n2,1)-x(n1,1);
      vv = .125*(dx*Ax+dy*Ay);
      for (int k=0; k<nq; k++){
	ds0     =              s  (n1,k)  +s  (n2,k);
	ds1     =-.5  *(dx   *(sxx(n1,0,k)+sxx(n2,0,k))   
		       +dy   *(sxx(n1,1,k)+sxx(n2,1,k)));
	ds2     =-.125*(dx*dx*(sxx(n1,2,k)+sxx(n2,2,k))
		       +dy*dy*(sxx(n1,4,k)+sxx(n2,4,k))
		    +2.*dx*dy*(sxx(n1,3,k)+sxx(n2,3,k)));
	r(n1,k) -= vv*(ds0+ds1+ds2);
	r(n2,k) -= vv*(ds0-ds1+ds2);
      }}
    sxx.deallocate();
  }
}
