#include "StrandBlockSolver.h"


void StrandBlockSolver::lhsViscousFine()
{
  int jj=nPstr+2,nn=nFaces+nBedges;
  Array4D<double> aa(nq,nq,jj,nn);
  for (int n=0; n<nFaces+nBedges; n++)
    for (int j=0; j<nPstr+2; j++)
      for (int k=0; k<nq; k++)
	for (int l=0; l<nq; l++) aa(l,k,j,n) = 0.;


  // obtain the Jacobian of the primitive variables
  int npts=jj*nn;
  Array4D<double> pj(nq,nq,jj,nn);
  sys->lhsPrimVarJacobian(npts,&q(0,0,0),&qa(0,0,0),&pj(0,0,0,0));


  // unstructured faces
  int c1,c2,n1,n2,fc,jm,jp,m,mm,m1l,m1u,m2l,m2u;
  double dx1,dy1,dx2,dy2,ds,a[nq*nq],b1[nq*nq],b2[nq*nq],B[nq],qe[nq],qae[nqa];
  npts = 1;
  for (int n=0; n<nEdges; n++){
    c1              = edge(0,n);
    c2              = edge(1,n);
    n1              = edgn(n);
    m1l             = ncsc(c1  );
    m1u             = ncsc(c1+1);
    m2l             = ncsc(c2  );
    m2u             = ncsc(c2+1);
    fc              = fClip(c1);
    if (fClip(c2) > fc) fc = fClip(c2);
  for (int j=1; j<fc+1; j++){
    jm              = j-1;
    dx1             = x (0,j,n1)-x (0,jm,n1);
    dy1             = x (1,j,n1)-x (1,jm,n1);
    dx2             = xc(0,j,c2)-xc(0,j ,c1);
    dy2             = xc(1,j,c2)-xc(1,j ,c1);
    ds              = 1./(dx1*dy2-dx2*dy1);
    B[0]            =-ds*dy1;
    B[1]            = ds*dx1;
    for (int k=0; k<nq ; k++) qe [k] = .5*(qp (k,j,n1)+qp (k,jm,n1));
    for (int k=0; k<nqa; k++) qae[k] = .5*(qap(k,j,n1)+qap(k,jm,n1));
    sys->lhsVisFluxJacobian(npts,&facs(0,j,n),&B[0],&qe[0],&qae[0],&a[0]);
    matmul(nq,nq,nq,&a[0],&pj(0,0,j,c1),&b1[0]);
    matmul(nq,nq,nq,&a[0],&pj(0,0,j,c2),&b2[0]);

    // diagonal contributions
    for (int k=0; k<nq; k++){
      m             = k*nq;
    for (int l=0; l<nq; l++){
      dd(l,k,j,c1) += b1[m+l];
      dd(l,k,j,c2) += b2[m+l];
    }}

    // off-diagonal contributions
    for (int k=0; k<nq; k++){
      m             = k*nq;
    for (int l=0; l<nq; l++){
      aa(l,k,j,c1)  = b1[m+l];
      aa(l,k,j,c2)  = b2[m+l];
    }}

    for (int mm=m1l; mm<m1u; mm++){
      nn            = csc(mm);
      for (int k=0; k<nq; k++){
	m           = k*nq;
      for (int l=0; l<nq; l++)
	bu(l,k,j,mm) -= aa(l,k,j,nn);
      }}

    for (int mm=m2l; mm<m2u; mm++){
      nn            = csc(mm);
      for (int k=0; k<nq; k++){
	m           = k*nq;
      for (int l=0; l<nq; l++)
	bu(l,k,j,mm) -= aa(l,k,j,nn);
      }}

    // reset working array to zero
    for (int k=0; k<nq; k++){
      m             = k*nq;
    for (int l=0; l<nq; l++){
      aa(l,k,j,c1)  = 0.;
      aa(l,k,j,c2)  = 0.;
    }}
  }}

  aa.deallocate();


  // structured faces
  double eps=1.e-14;
  npts = 1;
  for (int n=0; n<nFaces-nGfaces; n++){
    n1             = face(0,n);
    n2             = face(1,n);
  for (int j=0; j<fClip(n)+1; j++){
    jp             = j+1;
    dx1            = x (0,j ,n2)-x (0,j,n1);
    dy1            = x (1,j ,n2)-x (1,j,n1);
    dx2            = xc(0,jp,n )-xc(0,j,n );
    dy2            = xc(1,jp,n )-xc(1,j,n );
    ds             = dx1*dy2-dx2*dy1;
  if (fabs(ds) < eps)
    for (int k=0; k<nq*nq; k++){
      b1[k] = 0.;
      b2[k] = 0.;
    }
  else{
    ds             = 1./ds;
    B[0]           =-ds*dy1;
    B[1]           = ds*dx1;
    for (int k=0; k<nq;  k++) qe [k] = .5*(qp (k,j,n1)+qp (k,j,n2));
    for (int k=0; k<nqa; k++) qae[k] = .5*(qap(k,j,n1)+qap(k,j,n2));
    sys->lhsVisFluxJacobian(npts,&facu(0,j,n),&B[0],&qe[0],&qae[0],&a[0]);
    matmul(nq,nq,nq,&a[0],&pj(0,0,j ,n),&b1[0]);
    matmul(nq,nq,nq,&a[0],&pj(0,0,jp,n),&b2[0]);
  }
   for (int k=0; k<nq; k++){
     m             = k*nq;
   for (int l=0; l<nq; l++){
     dd(l,k,j ,n) += b1[m+l];
     dp(l,k,j ,n) -= b2[m+l];
     dd(l,k,jp,n) += b2[m+l];
     dm(l,k,jp,n) -= b1[m+l];
   }}}}

  pj.deallocate();
}
