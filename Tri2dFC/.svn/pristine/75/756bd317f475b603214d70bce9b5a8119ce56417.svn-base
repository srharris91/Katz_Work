#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::rhsViscous()
{
  if (order == 1 || order == 2){
    // Galerkin
    int n1,n2,n3,nn;
    double xa,ya,xb,yb,xc,yc,vr,fv[nq],gv[nq];
    Array2D<double> qc(3,nq),qac(3,nqa),cx(3,2);
    for (int n=0; n<nTri; n++){
      n1 = tri(n,0);
      n2 = tri(n,1);
      n3 = tri(n,2);
      xa = x(n1,0);
      ya = x(n1,1);
      xb = x(n2,0);
      yb = x(n2,1);
      xc = x(n3,0);
      yc = x(n3,1);
      vr = xa*(yb-yc)+xb*(yc-ya)+xc*(ya-yb);
      vr = 2./vr;
      
      for (int i=0; i<3; i++){
	for (int k=0; k<nq ; k++) qc (i,k) = q (n1,k);
	for (int k=0; k<nqa; k++) qac(i,k) = qa(n1,k);
	cx(i,0) = .5*(x(n3,1)-x(n2,1));
	cx(i,1) = .5*(x(n2,0)-x(n3,0));
	nn      = n1;
	n1      = n2;
	n2      = n3;
	n3      = nn;
      }
      
      sys->rhsVisFluxGalerkin(1,&vr,&cx(0,0),&qc(0,0),&qac(0,0),&fv[0],&gv[0]);
      for (int k=0; k<nq ; k++) d(n1,k) -=(cx(0,0)*fv[k]+cx(0,1)*gv[k]);
      for (int k=0; k<nq ; k++) d(n2,k) -=(cx(1,0)*fv[k]+cx(1,1)*gv[k]);
      for (int k=0; k<nq ; k++) d(n3,k) -=(cx(2,0)*fv[k]+cx(2,1)*gv[k]);
    }
    qc.deallocate();
    qac.deallocate();
    cx.deallocate();
  }


  else{// third-order scheme
    int m,l1,l2,n1,n2;
    double dx,dy,fk,Ax,Ay,cf,cg,a=.5;
    Array2D<double> qaxE(nne,nqaGradQa),qayE(nne,nqaGradQa);
    Array2D<double> f(nne,nq),g(nne,nq);
    Array2D<double> fx(nne,nq),gx(nne,nq),fy(nne,nq),gy(nne,nq);
    for (int n=0; n<nElem; n++){

      qaxE.set(0.);
      qayE.set(0.);
      fx.set(0.);
      fy.set(0.);
      gx.set(0.);
      gy.set(0.);

      // compute local gradients at the nodes in this element
      for (int i=0; i<nne; i++)
	for (int j=0; j<nne; j++){
	  dx = dxg(n,i,j,0);
	  dy = dxg(n,i,j,1);
	  m  = elem(n,j);
	  for (int k=0; k<nqaGradQa; k++){
	    qaxE(i,k) += qa(m,iqagrad(k))*dx;
	    qayE(i,k) += qa(m,iqagrad(k))*dy;
	  }}

      // compute viscous fluxes at the nodes in this element
      for (int i=0; i<nne; i++){
	m = elem(n,i);
	sys->rhsVisFlux(1,&q(m,0),&qa(m,0),&qaxE(i,0),&qayE(i,0),
			&f(i,0),&g(i,0));
      }

      // compute gradients of viscous fluxes
      for (int i=0; i<nne; i++)
	for (int j=0; j<nne; j++){
	  dx = dxg(n,i,j,0);
	  dy = dxg(n,i,j,1);
	  for (int k=0; k<nq ; k++){
	    fx(i,k) += f(j,k)*dx;
	    gx(i,k) += g(j,k)*dx;
	    fy(i,k) += f(j,k)*dy;
	    gy(i,k) += g(j,k)*dy;
	  }}

      // compute directed viscous fluxes and corrections
      // at edges and distribute to nodes
      for (int i=0; i<nee; i++){
	l1 = edgeE(i,0);
	l2 = edgeE(i,1);
	n1 = elem(n,l1);
	n2 = elem(n,l2);
	Ax = a*areaE(n,i,0);
	Ay = a*areaE(n,i,1);
	dx = .5*(x(n2,0)-x(n1,0));
	dy = .5*(x(n2,1)-x(n1,1));
	for (int k=0; k<nq; k++){
	  cf       = dx*(fx(l2,k)-fx(l1,k))+
		     dy*(fy(l2,k)-fy(l1,k));
	  cg       = dx*(gx(l2,k)-gx(l1,k))+
		     dy*(gy(l2,k)-gy(l1,k));
	  fk       = Ax*(f(l1,k)+f(l2,k)-cf)+Ay*(g(l1,k)+g(l2,k)-cg);
	  d(n1,k) -= fk;
	  d(n2,k) += fk;
	}}}

    qaxE.deallocate();
    qayE.deallocate();
    f.deallocate();
    g.deallocate();
    fx.deallocate();
    gx.deallocate();
    fy.deallocate();
    gy.deallocate();
  }
}
