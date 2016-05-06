#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::rhsTime()
{
  // form unsteady time derivative source term
  Array2D<double> t(nNode,nq);
  for (int n=0; n<nNode; n++)
    for (int k=0; k<nq; k++) t(n,k) = bdf(0)*q(n,k);
  for (int m=0; m<timeAcc; m++)
    for (int n=0; n<nNode; n++)
      for (int k=0; k<nq; k++) t(n,k) += bdf(m+1)*qt(m,n,k);
  for (int n=0; n<nNode; n++)
    for (int k=0; k<nq; k++) t(n,k) /= dtUnsteady;
  

  // discretize the source term based on the order of accuracy
  if      (order == 1)
    for (int n=0; n<nNode; n++)
      for (int k=0; k<nq; k++) r(n,k) += v(n)*t(n,k);

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
	a        = vv*(t(n1,k)+t(n2,k));
	r(n1,k) += a;
	r(n2,k) += a;
      }}}

  else if (order == 3){
    int n1,n2;
    double Ax,Ay,dx,dy,vv,dt0,dt1,dt2;
    Array3D<double> txx(nNode,5,nq);

    // compute gradient and Hessian of time derivative
    hessian(nq,&t(0,0),&txx(0,0,0));

    for (int n=0; n<nEdge; n++){
      n1 = edge(n,0);
      n2 = edge(n,1);
      Ax = area(n,0);
      Ay = area(n,1);
      dx = x(n2,0)-x(n1,0);
      dy = x(n2,1)-x(n1,1);
      vv = .25*(dx*Ax+dy*Ay);
      for (int k=0; k<nq; k++){
	dt0     =              t  (n1,k)  +t  (n2,k);
	dt1     =-.5  *(dx   *(txx(n1,0,k)+txx(n2,0,k))   
		       +dy   *(txx(n1,1,k)+txx(n2,1,k)));
	dt2     =-.125*(dx*dx*(txx(n1,2,k)+txx(n2,2,k))
		       +dy*dy*(txx(n1,4,k)+txx(n2,4,k))
		    +2.*dx*dy*(txx(n1,3,k)+txx(n2,3,k)));
	r(n1,k) += .5*(dt0+dt1+dt2)*vv;
	r(n2,k) += .5*(dt0-dt1+dt2)*vv;
      }}
    txx.deallocate();
  }

  t.deallocate();
}
