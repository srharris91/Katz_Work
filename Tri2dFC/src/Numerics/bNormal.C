#include "Tri2dFCBlockSolver.h"
#include "lagrangePoly.h"


void Tri2dFCBlockSolver::bNormal()
{
  // for each cell, compute Jacobian terms using the highest order
  // approximation available in the mesh
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
    }}


  // specify the local node numbers along each edge
  Array2D<int> keyE(3,orderE+1);
  if      (orderE == 1){
    keyE(0,0) = 0;
    keyE(0,1) = 1;
    keyE(1,0) = 1;
    keyE(1,1) = 2;
    keyE(2,0) = 2;
    keyE(2,1) = 0;
  }
  else if (orderE == 2){
    keyE(0,0) = 0;
    keyE(0,1) = 3;
    keyE(0,2) = 1;
    keyE(1,0) = 1;
    keyE(1,1) = 4;
    keyE(1,2) = 2;
    keyE(2,0) = 2;
    keyE(2,1) = 5;
    keyE(2,2) = 0;
  }
  else if (orderE == 3){
    keyE(0,0) = 0;
    keyE(0,1) = 3;
    keyE(0,2) = 4;
    keyE(0,3) = 1;
    keyE(1,0) = 1;
    keyE(1,1) = 5;
    keyE(1,2) = 6;
    keyE(1,3) = 2;
    keyE(2,0) = 2;
    keyE(2,1) = 7;
    keyE(2,2) = 8;
    keyE(2,3) = 0;
  }


  // local normals for each element edge
  if (orderE == 1) nsq = 1; // number of surface quadrature points
  else nsq = 2;
  Array2D<double> gp(nsq,2); // gauss points and weights [-1,1]
  if (nsq == 1){
    gp(0,0) = 0.;
    gp(0,1) = 2.;
  }
  else{
    gp(0,0) =-1./sqrt(3.);
    gp(0,1) = 1.;
    gp(1,0) = 1./sqrt(3.);
    gp(1,1) = 1.;
  }
  Array2D<double> sn(3,2);
  Array3D<double> sq(3,nsq,4);
  double r0,s0,r1,s1,a;
  int n1,n2,m;
  j       = 0;
  r0      =-1.;
  s0      =-1./sqrt(3.);
  r1      = 1.;
  s1      =-1./sqrt(3.);
  sn(j,0) = s1-s0;
  sn(j,1) = r0-r1;
  a       = 1./sqrt(pow(sn(j,0),2)+pow(sn(j,1),2));
  sn(j,0) = sn(j,0)*a;
  sn(j,1) = sn(j,1)*a;
  for (int n=0; n<nsq; n++){
    sq(j,n,0) = r0+.5*(gp(n,0)+1.)*(r1-r0); //r quadrature loc.
    sq(j,n,1) = s0+.5*(gp(n,1)+1.)*(s1-s0); //s quadrature loc.
    sq(j,n,2) = 1.*gp(n,2); //w_j*e
    sq(j,n,3) =-1./sqrt(3.); //b
  }

  j       = 1;
  r0      = r1;
  s0      = s1;
  r1      = 0.;
  s1      = 2./sqrt(3.);
  sn(j,0) = s1-s0;
  sn(j,1) = r0-r1;
  a       = 1./sqrt(pow(sn(j,0),2)+pow(sn(j,1),2));
  sn(j,0) = sn(j,0)*a;
  sn(j,1) = sn(j,1)*a;
  for (int n=0; n<nsq; n++){
    sq(j,n,0) = r0+.5*(gp(n,0)+1.)*(r1-r0);
    sq(j,n,1) = s0+.5*(gp(n,1)+1.)*(s1-s0);
    sq(j,n,2) =-.5*gp(n,2); //w_j*e
    sq(j,n,3) =-3./sqrt(3.); //b
  }

  j       = 2;
  r0      = r1;
  s0      = s1;
  r1      =-1.;
  s1      =-1./sqrt(3.);
  sn(j,0) = s1-s0;
  sn(j,1) = r0-r1;
  a       = 1./sqrt(pow(sn(j,0),2)+pow(sn(j,1),2));
  sn(j,0) = sn(j,0)*a;
  sn(j,1) = sn(j,1)*a;
  for (int n=0; n<nsq; n++){
    sq(j,n,0) = r0+.5*(gp(n,0)+1.)*(r1-r0);
    sq(j,n,1) = s0+.5*(gp(n,1)+1.)*(s1-s0);
    sq(j,n,2) =-.5*gp(n,2); //w_j*e
    sq(j,n,3) = 3./sqrt(3.); //b
  }


  // find the boundary tag for each element edge
  // flag boundary nodes with tags of two boundary edge numbers that share it
  Array2D<int> nflag(nNode,2);
  Array2D<int> eTag(nElem,3);
  nflag.set(-1);
  for (int n=nEdge-nEdgeBd; n<nEdge; n++){
    n1 = edge(n,0);
    n2 = edge(n,1);
    if (nflag(n1,0) == -1) nflag(n1,0) = edge(n,3);
    else nflag(n1,1) = edge(n,3);
    if (nflag(n2,0) == -1) nflag(n2,0) = edge(n,3);
    else nflag(n2,1) = edge(n,3);
  }

  int k1;
  for (int n=0; n<nElem; n++)
    for (int j=0; j<3; j++){
      n1 = elem(n,keyE(j,0     )); // first node on element edge j
      n2 = elem(n,keyE(j,orderE)); // last  node on element edge j
      if (nflag(n1,0) == -1 || nflag(n2,0) == -1) k1 =-1;
      else{
	if      (nflag(n1,0) == nflag(n2,0) ||
	         nflag(n1,0) == nflag(n2,1)) k1 = nflag(n1,0);
	else if (nflag(n1,1) == nflag(n2,0) ||
		 nflag(n1,1) == nflag(n2,1)) k1 = nflag(n1,1);
	else k1 =-1;
      }
      eTag(n,j) = k1;
      //cout << n << " " << n1 << " " << n2 << " " << k1 << endl;
    }

  ln.set(0.);
  int k;
  double xrn,yrn,xsn,ysn,jcn,dx,dy;
  for (int n=0; n<nElem; n++)
    for (int j=0; j<3; j++)
      if (eTag(n,j) >= 0)
	for (int i=0; i<orderE+1; i++){
	  k = elem(n,keyE(j,i));
	  m = k-nNode+nNodeBd;
	  if (nodeBd(m) == eTag(n,j)){
	    xrn      = xr(n,keyE(j,i));
	    yrn      = yr(n,keyE(j,i));
	    xsn      = xs(n,keyE(j,i));
	    ysn      = ys(n,keyE(j,i));
	    dx       = xrn*sn(j,1)-xsn*sn(j,0);
	    dy       = yrn*sn(j,1)-ysn*sn(j,0);
	    a        = 1./sqrt(dx*dx+dy*dy);
	    ln(m,0) +=-dy*a;
	    ln(m,1) += dx*a;
	  }}

  for (int n=0; n<nNodeBd; n++){
    a        = 1./sqrt(ln(n,0)*ln(n,0)+ln(n,1)*ln(n,1));
    ln(n,0) *= a;
    ln(n,1) *= a;
  }

  /*
  m = 0;
  for (int n=nNode-nNodeBd; n<nNode; n++){
    cout << n << " " << m << " "
         << nodeBd(m) << " " << ln(m,0) << " " << ln(m,1) << endl;
    m++;
  }
  */


  // detect sharp corners and set them to the "nothing" condition
  if (sourceMMS == 0){
    int n1,n2,p1,p2;
    double dx,dy,ds,nx,ny;
    Array2D<double> nrm(nNodeBd,4);
    nrm.set(0.);
    for (int n=nEdge-nEdgeBd; n<nEdge; n++){
      n1 = edge(n,0);
      n2 = edge(n,1);
      dx = x(n2,1)-x(n1,1);
      dy = x(n1,0)-x(n2,0);
      ds = 1./sqrt(dx*dx+dy*dy);
      dx = dx*ds;
      dy = dy*ds;
      p1 = n1-nNode+nNodeBd;
      p2 = n2-nNode+nNodeBd;
      nx = nrm(p1,0);
      ny = nrm(p1,1);
      ds = nx*nx+ny*ny;
      if (ds < .5){
	nrm(p1,0) = dx;
	nrm(p1,1) = dy;
      }
      else{
	nrm(p1,2) = dx;
	nrm(p1,3) = dy;
      }
      nx = nrm(p2,0);
      ny = nrm(p2,1);
      ds = nx*nx+ny*ny;
      if (ds < .5){
	nrm(p2,0) = dx;
	nrm(p2,1) = dy;
      }
      else{
	nrm(p2,2) = dx;
	nrm(p2,3) = dy;
      }
    }

    double pi=4.*atan(1.),t=91./180.*pi; //sharp corner at 91 deg.
    for (int n=0; n<nNodeBd; n++){
      dx = nrm(n,0);
      dy = nrm(n,1);
      nx = nrm(n,2);
      ny = nrm(n,3);
      ds = dx*nx+dy*ny;
      if (ds < cos(t)) nodeBd(n) = nCompBd-1;
    }
    nrm.deallocate();
  }


  // form element edges and data structures for surface force integration
  nEdgeQ = 0;
  for (int n=0; n<nElem; n++)
    for (int j=0; j<3; j++) 
      if (eTag(n,j) != -1) nEdgeQ++;
  edgeQ.allocate(nEdgeQ,3);
  nxQ.allocate(nEdgeQ,nsq,2);
  gxS.allocate(nEdgeQ,nsq,nne,3);

  Array2D<double> xrQ(nEdgeQ,nsq),yrQ(nEdgeQ,nsq),xsQ(nEdgeQ,nsq),
    ysQ(nEdgeQ,nsq),jacQ(nEdgeQ,nsq);
  Array3D<double> lQ(3,nsq,nne),lrQ(3,nsq,nne),lsQ(3,nsq,nne);
  xrQ.set(0.);
  yrQ.set(0.);
  xsQ.set(0.);
  ysQ.set(0.);
  lQ.set(0.);
  lrQ.set(0.);
  lsQ.set(0.);

  // lrQ(i,j) = (dl_j/dr)_i (a row is all Lagrange polynomials (derivatives)
  // evaluated at a single surface quadrature point i,
  // same with the other derivatives)
  int jj;
  for (int n=0; n<nne; n++) // nth Lagrange polynomial
    for (int j=0; j<3; j++){
      for (int i=0; i<nsq; i++){ // ith surface quadrature point
	jj = 0;
	ri = sq(j,i,0);
	si = sq(j,i,1);
	for (int k=0; k<=orderE; k++)
	  for (int l=0; l<=orderE-k; l++){
	    km          = max(0,k-1);
	    lm          = max(0,l-1);
	    lQ (j,i,n) +=            pow(ri,k )*pow(si,l )*lc(n,jj  );
	    lrQ(j,i,n) +=((double)k)*pow(ri,km)*pow(si,l )*lc(n,jj  );
	    lsQ(j,i,n) +=((double)l)*pow(ri,k )*pow(si,lm)*lc(n,jj++);
	  }}}

  // Jacobian terms
  m = 0;
  double wi,bi;
  for (int n=0; n<nElem; n++)
    for (int j=0; j<3; j++) 
      if (eTag(n,j) != -1){
	edgeQ(m,0) = n;
	edgeQ(m,1) = j;
	edgeQ(m,2) = eTag(n,j);
	for (int i=0; i<nsq; i++){
	  ri = sq(j,i,0); // location and weight on the edge
	  si = sq(j,i,1);
	  wi = sq(j,i,2);
	  bi = sq(j,i,3);
	  for (int mm=0; mm<nne; mm++){ //mmth Lagrange polynomial
	    jj  = 0;
	    lrm = 0.;
	    lsm = 0.;
	    for (int k=0; k<=orderE; k++){
	      for (int l=0; l<=orderE-k; l++){
		km   = max(0,k-1);
		lm   = max(0,l-1);
		lrm +=((double)k)*pow(ri,km)*pow(si,l )*lc(mm,jj  );
		lsm +=((double)l)*pow(ri,k )*pow(si,lm)*lc(mm,jj++);
	      }}
	    xrQ(m,i) += lrm*x(elem(n,mm),0);
	    yrQ(m,i) += lrm*x(elem(n,mm),1);
	    xsQ(m,i) += lsm*x(elem(n,mm),0);
	    ysQ(m,i) += lsm*x(elem(n,mm),1);
	  }
	  jacQ(m,i)  = xrQ(m,i)*ysQ(m,i)-yrQ(m,i)*xsQ(m,i);
	  nxQ(m,i,0) = wi*(yrQ(m,i)+bi*ysQ(m,i));
	  nxQ(m,i,1) =-wi*(xrQ(m,i)+bi*xsQ(m,i));
	}
	m++;
      }

  /*
  for (int n=0; n<nEdgeQ; n++)
    cout << n << " " << edgeQ(n,0) << " " << edgeQ(n,1) << endl;
  exit(0);
  */

  double jci,xri,yri,xsi,ysi;
  for (int n=0; n<nEdgeQ; n++){
    jj  = edgeQ(n,1); // edge number on the element
    for (int i=0; i<nsq; i++){ //for each surface quadrature point
      jci = 1./jacQ(n,i);
      xri = xrQ(n,i)*jci;
      yri = yrQ(n,i)*jci;
      xsi = xsQ(n,i)*jci;
      ysi = ysQ(n,i)*jci;
      for (int j=0; j<nne; j++){
	gxS(n,i,j,0) = lQ (jj,i,j);
	gxS(n,i,j,1) = lrQ(jj,i,j)*ysi-lsQ(jj,i,j)*yri;
	gxS(n,i,j,2) =-lrQ(jj,i,j)*xsi+lsQ(jj,i,j)*xri;
      }}}
  

  // deallocate work arrays
  xr.deallocate();
  yr.deallocate();
  xs.deallocate();
  ys.deallocate();
  jac.deallocate();
  rs.deallocate();
  lc.deallocate();
  keyE.deallocate();
  gp.deallocate();
  sn.deallocate();
  sq.deallocate();
  nflag.deallocate();
  eTag.deallocate();
  xrQ.deallocate();
  yrQ.deallocate();
  xsQ.deallocate();
  ysQ.deallocate();
  jacQ.deallocate();
  lQ.deallocate();
  lrQ.deallocate();
  lsQ.deallocate();
  lc.deallocate();
}
