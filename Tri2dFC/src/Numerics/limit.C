#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::limit()
{
  if (limiter == 0) lim.set(1.);

  else{
    int n1,n2,n3,t1,t2;
    double xa,ya,xb,yb,xc,yc,dxa,dya,dxb,dyb,dxc,dyc,vr,dx,dy,dql,dqr,a;
    Array3D<double> qxt;
    qxt.allocate(nTri,nq,2);
    for (int n=0; n<nTri; n++){
      n1  = tri(n,0);
      n2  = tri(n,1);
      n3  = tri(n,2);
      xa  = x(n1,0);
      ya  = x(n1,1);
      xb  = x(n2,0);
      yb  = x(n2,1);
      xc  = x(n3,0);
      yc  = x(n3,1);
      dxa = xb-xc;
      dxb = xc-xa;
      dxc = xa-xb;
      dya = yb-yc;
      dyb = yc-ya;
      dyc = ya-yb;
      vr  = 1./(xa*dya+xb*dyb+xc*dyc);
      for (int k=0; k<nq; k++){
	qxt(n,k,0) = (q(n1,k)*dya+q(n2,k)*dyb+q(n3,k)*dyc)*vr;
	qxt(n,k,1) =-(q(n1,k)*dxa+q(n2,k)*dxb+q(n3,k)*dxc)*vr;
      }}

    for (int n=0; n<nEdge; n++){
      n1 = edge(n,0);
      n2 = edge(n,1);
      t1 = edge(n,4);
      t2 = edge(n,5);
      dx = x(n2,0)-x(n1,0);
      dy = x(n2,1)-x(n1,1);
      for (int k=0; k<nq; k++){
	dql      = dx*qxt(t1,k,0)+dy*qxt(t1,k,1);
	dqr      = dx*qxt(t2,k,0)+dy*qxt(t2,k,1);
	a        =(dqr-dql)/max(fabs(dqr)+fabs(dql),dlim(k));
	lim(n,k) = 1.-pow(fabs(a),3);
      }}

    /*
    double minlim=1.;
    for (int n=0; n<nEdge; n++) if (lim(n,0) < minlim) minlim = lim(n,0);
    cout << minlim << endl;
    */

    qxt.deallocate();
  }
}
