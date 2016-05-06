#include "StrandBlockSolver.h"


void StrandBlockSolver::lhsViscousCoarse()
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
  int c1,c2,fc,jp,m,mm,m1l,m1u,m2l,m2u;
  double ve,a[nq*nq],b1[nq*nq],b2[nq*nq],qe[nq],qae[nqa];
  npts = 1;
  for (int n=0; n<nEdges; n++){
    c1              = edge(0,n);
    c2              = edge(1,n);
    m1l             = ncsc(c1  );
    m1u             = ncsc(c1+1);
    m2l             = ncsc(c2  );
    m2u             = ncsc(c2+1);
    fc              = fClip(c1);
    if (fClip(c2) > fc) fc = fClip(c2);
  for (int j=1; j<fc+1; j++){
    for (int k=0; k<nq ; k++) qe [k] = .5*(q (k,j,c1)+q (k,j,c2));
    for (int k=0; k<nqa; k++) qae[k] = .5*(qa(k,j,c1)+qa(k,j,c2));
    ve = 2./(v(j,c1)+v(j,c2));
    sys->lhsVisFluxJacobian(npts,&facs(0,j,n),&facs(0,j,n),
			    &qe[0],&qae[0],&a[0]);
    for (int k=0; k<nq*nq; k++) a[k] *= ve;
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
  npts = 1;
  for (int n=0; n<nFaces-nGfaces; n++){
  for (int j=0; j<fClip(n)+1; j++){
    jp             = j+1;
    for (int k=0; k<nq;  k++) qe [k] = .5*(q (k,j,n)+q (k,jp,n));
    for (int k=0; k<nqa; k++) qae[k] = .5*(qa(k,j,n)+qa(k,jp,n));
    ve = 2./(v(j,n)+v(jp,n));
    sys->lhsVisFluxJacobian(npts,&facu(0,j,n),&facu(0,j,n),
			    &qe[0],&qae[0],&a[0]);
    matmul(nq,nq,nq,&a[0],&pj(0,0,j ,n),&b1[0]);
    matmul(nq,nq,nq,&a[0],&pj(0,0,jp,n),&b2[0]);
   for (int k=0; k<nq; k++){
     m             = k*nq;
   for (int l=0; l<nq; l++){
     dd(l,k,j ,n) += b1[m+l]*ve;
     dp(l,k,j ,n) -= b2[m+l]*ve;
     dd(l,k,jp,n) += b2[m+l]*ve;
     dm(l,k,jp,n) -= b1[m+l]*ve;
   }}}}

  pj.deallocate();
}
