#include "Tri2dFCBlockSolver.h"
#include "lagrangePoly.h"


void Tri2dFCBlockSolver::gradSetupCubic()
{
  // for each cell, compute Jacobian terms using the highest order
  // approximation available in the global mesh
  Array2D <double> xr(nElem,nne),yr(nElem,nne),xs(nElem,nne),
    ys(nElem,nne),jac(nElem,nne),rs(nne,3),lc(nne,nne);
  xr.set(0.);
  yr.set(0.);
  xs.set(0.);
  ys.set(0.);
  jac.set(0.);

  int orderE=3-level;
  solutionPoints(orderE,
		 spacing,
		 &rs(0,0));

  bool test=true;
  lagrangePoly(test,
	       orderE,
	       &rs(0,0),
	       &lc(0,0));

  int j,km,lm;
  double dL,L0,L1,lrm,lsm,ri,si;
  for (int n=0; n<nElem; n++){

    // evaluate the Jacobian terms at the mesh points
    for (int i=0; i<nne; i++){ //ith mesh point
      ri = rs(i,0);
      si = rs(i,1);
      for (int m=0; m<nne; m++){ //mth Lagrange polynomial used in the mapping
	j   = 0;
	lrm = 0.;
	lsm = 0.;
	for (int k=0; k<=orderE; k++)
	  for (int l=0; l<=orderE-k; l++){
	    km   = max(0,k-1);
	    lm   = max(0,l-1);
	    lrm +=((double)k)*pow(ri,km)*pow(si,l )*lc(m,j  );
	    lsm +=((double)l)*pow(ri,k )*pow(si,lm)*lc(m,j++);
	  }
	xr(n,i) += lrm*x(elem(n,m),0);
	yr(n,i) += lrm*x(elem(n,m),1);
	xs(n,i) += lsm*x(elem(n,m),0);
	ys(n,i) += lsm*x(elem(n,m),1);
      }
      jac(n,i) = xr(n,i)*ys(n,i)-yr(n,i)*xs(n,i);
      //cout << n << " " << elem(n,i) << " " << jac(n,i) << endl;
    }}
  //exit(0);

  // lr(i,j) = (dl_j/dr)_i (a row is all Lagrange polynomials (derivatives)
  // evaluated at a single mesh point i, same with the other derivatives)
  Array2D<double> lr(nne,nne),ls(nne,nne),lrr(nne,nne),lss(nne,nne),
    lrs(nne,nne);
  lr.set(0.);
  ls.set(0.);
  lrr.set(0.);
  lss.set(0.);
  lrs.set(0.);
  int kmm,lmm;
  for (int n=0; n<nne; n++) // nth Lagrange polynomial
    for (int i=0; i<nne; i++){ // ith mesh point
      j  = 0;
      ri = rs(i,0);
      si = rs(i,1);
      for (int k=0; k<=orderE; k++)
	for (int l=0; l<=orderE-k; l++){
	  km        = max(0,k-1);
	  lm        = max(0,l-1);
	  kmm       = max(0,k-2);
	  lmm       = max(0,l-2);
	  lr (i,n) +=((double)k)*pow(ri,km)*pow(si,l )*lc(n,j);
	  ls (i,n) +=((double)l)*pow(ri,k )*pow(si,lm)*lc(n,j);
	  lrr(i,n) +=((double)(k*km))*pow(ri,kmm)*pow(si,l  )*lc(n,j  );
	  lss(i,n) +=((double)(l*lm))*pow(ri,k  )*pow(si,lmm)*lc(n,j  );
	  lrs(i,n) +=((double)(k*l ))*pow(ri,km )*pow(si,lm )*lc(n,j++);
	}}


  // cubic gradient coefficients for the viscous terms
  double jci,xri,yri,xsi,ysi,dx,dy;
  for (int n=0; n<nElem; n++)
    for (int i=0; i<nne; i++){
      jci = jac(n,i);
      xri = xr(n,i);
      yri = yr(n,i);
      xsi = xs(n,i);
      ysi = ys(n,i);
      for (int j=0; j<nne; j++){
	dxg(n,i,j,0) =( lr(i,j)*ysi-ls(i,j)*yri)/jci;
	dxg(n,i,j,1) =(-lr(i,j)*xsi+ls(i,j)*xri)/jci;
      }}


  // compute fully cubic FEM gradient coefficients
  Array1D<double> sumj(nNode);
  sumj.set(0.);
  for (int n=0; n<nElem; n++)
    for (int i=0; i<nne; i++) sumj(elem(n,i)) += jac(n,i);
  for (int n=0; n<nNode; n++) sumj(n) = 1./sumj(n);

  int k1,k2;
  Array2D<double> ax(nNode,2);
  ax.set(0.);
  gxC.set(0.);
  for (int n=0; n<nElem; n++)
    for (int i=0; i<nne; i++){
      k1  = elem(n,i);
      xri = xr(n,i);
      yri = yr(n,i);
      xsi = xs(n,i);
      ysi = ys(n,i);
      for (int j=0; j<nne; j++){
	k2       = elem(n,j);
	ax(k2,0) = lr(i,j)*ysi-ls(i,j)*yri;
	ax(k2,1) =-lr(i,j)*xsi+ls(i,j)*xri;
      }
      for(int j=psp2(k1); j<psp2(k1+1); j++){
	k2        = psp1(j);
	gxC(j,0) += ax(k2,0);
	gxC(j,1) += ax(k2,1);
      }
      for (int j=0; j<nne; j++){
	k2       = elem(n,j);
	ax(k2,0) = 0.;
	ax(k2,1) = 0.;
      }}

  for (int n=0; n<nNode; n++)
    for (int i=psp2(n); i<psp2(n+1); i++){
      gxC(i,0) *= sumj(n);
      gxC(i,1) *= sumj(n);
    }


  // compute fully cubic FEM Hessian coefficients
  Array2D<double> axx(nNode,3);
  Array3D<double> px(nne,nne,2);
  axx.set(0.);
  gxx.set(0.);
  for (int n=0; n<nElem; n++){
    for (int i=0; i<nne; i++){
      jci = 1./jac(n,i);
      xri = xr(n,i)*jci;
      yri = yr(n,i)*jci;
      xsi = xs(n,i)*jci;
      ysi = ys(n,i)*jci;
      for (int j=0; j<nne; j++){
	px(i,j,0) = lr(i,j)*ysi-ls(i,j)*yri;
	px(i,j,1) =-lr(i,j)*xsi+ls(i,j)*xri;
      }}
    for (int i=0; i<nne; i++){
      k1  = elem(n,i);
      xri = xr(n,i);
      yri = yr(n,i);
      xsi = xs(n,i);
      ysi = ys(n,i);
      for (int j=0; j<nne; j++){
	dx = lr(i,j)*ysi-ls(i,j)*yri;
	dy =-lr(i,j)*xsi+ls(i,j)*xri;
	for (int k=0; k<nne; k++){
	  k2         = elem(n,k);
	  axx(k2,0) += px(j,k,0)*dx;
	  axx(k2,1) += px(j,k,0)*dy;
	  axx(k2,2) += px(j,k,1)*dy;
	}}
      for(int j=psp2(k1); j<psp2(k1+1); j++){
	k2        = psp1(j);
	gxx(j,0) += axx(k2,0);
	gxx(j,1) += axx(k2,1);
	gxx(j,2) += axx(k2,2);
      }
      for (int j=0; j<nne; j++){
	k2        = elem(n,j);
	axx(k2,0) = 0.;
	axx(k2,1) = 0.;
	axx(k2,2) = 0.;
      }}}

  for (int n=0; n<nNode; n++)
    for (int i=psp2(n); i<psp2(n+1); i++){
      gxx(i,0) *= sumj(n);
      gxx(i,1) *= sumj(n);
      gxx(i,2) *= sumj(n);
    }


  // deallocate work arrays
  xr.deallocate();
  yr.deallocate();
  xs.deallocate();
  ys.deallocate();
  jac.deallocate();
  rs.deallocate();
  lc.deallocate();
  lr.deallocate();
  ls.deallocate();
  lrr.deallocate();
  lss.deallocate();
  lrs.deallocate();
  sumj.deallocate();
  ax.deallocate();
  axx.deallocate();
  px.deallocate();
}
