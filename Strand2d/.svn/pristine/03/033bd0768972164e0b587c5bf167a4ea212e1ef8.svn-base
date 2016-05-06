
#include "StrandBlockSolver.h"


void StrandBlockSolver::lhsInviscid()
{
  int jj=nPstr+2,nn=nFaces+nBedges;
  Array4D<double> aa(nq,nq,jj,nn);
  for (int n=0; n<nFaces+nBedges; n++)
    for (int j=0; j<nPstr+2; j++)
      for (int k=0; k<nq; k++)
	for (int l=0; l<nq; l++) aa(l,k,j,n) = 0.;


  // unstructured faces
  int c1,c2,n1,n2,jm,jp,fc,m,mm,m1l,m1u,m2l,m2u,npts=1;
  double a1[nq*nq],a2[nq*nq],am[nq*nq],ap[nq*nq],Ax,Ay,ds,eps=1.e-14;
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
    sys->lhsInvFluxJacobian(npts,&facs(0,j,n),&xvs(j,n),
			    &q(0,j,c1),&qa(0,j,c1),&a1[0]);
    sys->lhsInvFluxJacobian(npts,&facs(0,j,n),&xvs(j,n),
			    &q(0,j,c2),&qa(0,j,c2),&a2[0]);

    for (int k=0; k<nq*nq; k++) a1[k] *= .5;
    for (int k=0; k<nq*nq; k++) a2[k] *= .5;

    // off-diagonal contributions
    for (int k=0; k<nq; k++){
      m             = k*nq;
    for (int l=0; l<nq; l++){
      aa(l,k,j,c1)  = a1[m+l];
      aa(l,k,j,c2)  = a2[m+l];
    }}

    for (int mm=m1l; mm<m1u; mm++){
      nn            = csc(mm);
      for (int k=0; k<nq; k++){
	m           = k*nq;
      for (int l=0; l<nq; l++)
	bu(l,k,j,mm) += aa(l,k,j,nn);
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
  for (int n=0; n<nFaces-nGfaces; n++){
  for (int j=1; j<fClip(n)+1; j++){
    jm           = j-1;
    jp           = j+1;
    Ax           = facu(0,jm,n);
    Ay           = facu(1,jm,n);
    ds           = Ax*Ax+Ay*Ay;
    if (fabs(ds) < eps){
      for (int k=0; k<nq; k++){
	m        = k*nq;
      for (int l=0; l<nq; l++)
	am[m+l]  = 0.;
      }}
    else
      sys->lhsInvFluxJacobian(npts,&facu(0,jm,n),&xvu(jm,n),
			      &q(0,jm,n),&qa(0,jm,n),&am[0]);
    sys->lhsInvFluxJacobian(npts,&facu(0,j ,n),&xvu(j ,n),
			    &q(0,jp,n),&qa(0,jp,n),&ap[0]);
    for (int k=0; k<nq; k++){
      m            = k*nq;
    for (int l=0; l<nq; l++){
      dm(l,k,j,n) -= (am[m+l]*.5);
      dp(l,k,j,n) += (ap[m+l]*.5);
    }}}}
}
