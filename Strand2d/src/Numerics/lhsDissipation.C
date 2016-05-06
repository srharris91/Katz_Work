#include "StrandBlockSolver.h"


void StrandBlockSolver::lhsDissipation()
{
  int jj=nPstr+2,nn=nFaces+nBedges;
  Array4D<double> aa(nq,nq,jj,nn);
  for (int n=0; n<nFaces+nBedges; n++)
    for (int j=0; j<nPstr+2; j++)
      for (int k=0; k<nq; k++)
	for (int l=0; l<nq; l++) aa(l,k,j,n) = 0.;


  // unstructured faces
  int c1,c2,n1,n2,jp,fc,m,mm,m1l,m1u,m2l,m2u,npts=1;
  double a[nq*nq];
  for (int n=0; n<nEdges; n++){
    c1            = edge(0,n);
    c2            = edge(1,n);
    m1l           = ncsc(c1  );
    m1u           = ncsc(c1+1);
    m2l           = ncsc(c2  );
    m2u           = ncsc(c2+1);
    fc            = fClip(c1);
    if (fClip(c2) > fc) fc = fClip(c2);
  for (int j=1; j<fc+1; j++){
    sys->lhsDisFluxJacobian(npts,&facs(0,j,n),&xvs(j,n),
			    &q(0,j,c1),&q(0,j,c2),&a[0]);

    for (int k=0; k<nq*nq; k++) a[k] *= .5;

    // diagonal contributions
    for (int k=0; k<nq; k++){
      m             = k*nq;
    for (int l=0; l<nq; l++){
      dd(l,k,j,c1) += a[m+l];
      dd(l,k,j,c2) += a[m+l];
    }}

    // off-diagonal contributions
    for (int k=0; k<nq; k++){
      m             = k*nq;
    for (int l=0; l<nq; l++){
      aa(l,k,j,c1)  = a[m+l];
      aa(l,k,j,c2)  = a[m+l];
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
  double Ax,Ay,ds,eps=1.e-14;
  for (int n=0; n<nFaces-nGfaces; n++){
    n1            = face(0,n);
    n2            = face(1,n);
  for (int j=0; j<fClip(n)+1; j++){
    jp            = j+1;
    Ax            = facu(0,j,n);
    Ay            = facu(1,j,n);
    ds            = Ax*Ax+Ay*Ay;
    if (fabs(ds) < eps){
      for (int k=0; k<nq; k++){
	m         = k*nq;
      for (int l=0; l<nq; l++)
	a[m+l]    = 0.;
      }}
    else
      sys->lhsDisFluxJacobian(npts,&facu(0,j,n),&xvu(j,n),
			      &q(0,j,n),&q(0,jp,n),&a[0]);
    for (int k=0; k<nq; k++){
      m             = k*nq;
    for (int l=0; l<nq; l++){
      dd(l,k,j ,n) += (a[m+l]*.5);
      dp(l,k,j ,n) -= (a[m+l]*.5);
      dd(l,k,jp,n) += (a[m+l]*.5);
      dm(l,k,jp,n) -= (a[m+l]*.5);
    }}}}
}
