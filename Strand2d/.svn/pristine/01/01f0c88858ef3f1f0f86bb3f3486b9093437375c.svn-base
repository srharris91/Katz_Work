#include "StrandBlockSolver.h"


void StrandBlockSolver::lhsBoundary(const int& mglevel)
{
  // initialize LHS linearization arrays to 0.,
  // setting corner and ghost points to the identity matrix
  int m1,m2;
  for (int n=0; n<nFaces+nBedges; n++){
    m1            = ncsc(n  );
    m2            = ncsc(n+1);
    for (int m=m1; m<m2; m++)
      for (int k=0; k<nq; k++)
      for (int l=0; l<nq; l++)
	bu(l,k,0,m) = 0.;

    for (int k=0; k<nq; k++)
    for (int l=0; l<nq; l++){
      dd(l,k,0,n) = 0.;
      dp(l,k,0,n) = 0.;
      dm(l,k,0,n) = 0.;
    }
    for (int j=fClip(n)+1; j<nPstr+2; j++){
      for (int k=0; k<nq; k++)
      for (int l=0; l<nq; l++){
	dd(l,k,j,n) = 0.;
	dp(l,k,j,n) = 0.;
	dm(l,k,j,n) = 0.;
      }
      for (int m=m1; m<m2; m++)
	for (int k=0; k<nq; k++)
	for (int l=0; l<nq; l++)
	  bu(l,k,j,m) = 0.;
    }
  }

  for (int n=nFaces-nGfaces; n<nFaces+nBedges; n++){
    m1            = ncsc(n  );
    m2            = ncsc(n+1);
    for (int m=m1; m<m2; m++){
      for (int j=0; j<nPstr+2; j++)
	for (int k=0; k<nq; k++)
	for (int l=0; l<nq; l++)
	  bu(l,k,j,m) = 0.;
    }
    for (int j=0; j<nPstr+2; j++){
      for (int k=0; k<nq; k++)
      for (int l=0; l<nq; l++){
	dd(l,k,j,n) = 0.;
	dp(l,k,j,n) = 0.;
	dm(l,k,j,n) = 0.;
      }}
  }

  for (int n=nFaces; n<nFaces+nBedges; n++){
    for (int k=0; k<nq; k++) dd(k,k,0,n) = 1.;
    for (int j=fClip(n)+1; j<nPstr+2; j++){
    for (int k=0; k<nq; k++) dd(k,k,j,n) = 1.;
    }}

  for (int n=nFaces-nGfaces; n<nFaces; n++){
    for (int j=0; j<nPstr+2; j++)
    for (int k=0; k<nq; k++) dd(k,k,j,n) = 1.;
  }


  int jj=nPstr+2,nn=nFaces+nBedges;
  Array4D<double> aa(nq,nq,jj,nn);
  for (int n=0; n<nFaces+nBedges; n++)
    for (int j=0; j<nPstr+2; j++)
      for (int k=0; k<nq; k++)
	for (int l=0; l<nq; l++) aa(l,k,j,n) = 0.;


  // unstructured boundaries
  int c1,c2,m,mb=-1,i,j,jp,npts=1;
  double nx[ndim],ds,a[nq*nq],ai[nq*nq];
  for (int n=nEdges-nBedges; n<nEdges; n++){
    mb++;
    c1           = edge(0,n);
    c2           = edge(1,n);
    m1           = ncsc(c2  );
    m2           = ncsc(c2+1);
  for (int j=1; j<fClip(c1)+1; j++){
    nx[0]        = facs(0,j,n);
    nx[1]        = facs(1,j,n);
    ds           = 1./sqrt(nx[0]*nx[0]+nx[1]*nx[1]);
    nx[0]       *= ds;
    nx[1]       *= ds;
    sys->lhsBCVectorSelfJacobian(npts,&bTag(mb),&nx[0],&q(0,j,c1),&qa(0,j,c1),
				 &q(0,j,c2),&qa(0,j,c2),&a[0]);
    sys->lhsBCVectorInteriorJacobian(npts,&bTag(mb),&nx[0],&q(0,j,c1),
				     &qa(0,j,c1),&q(0,j,c2),&qa(0,j,c2),&ai[0]);

    // diagonal contributions
    for (int k=0; k<nq; k++){
      i            = k*nq;
    for (int l=0; l<nq; l++)
      dd(l,k,j,c2) = a[i+l];
    }

    // off-diagonal contributions
    for (int k=0; k<nq; k++){
      m             = k*nq;
    for (int l=0; l<nq; l++)
      aa(l,k,j,c1)  = ai[m+l];
    }

    for (int mm=m1; mm<m2; mm++){
      nn            = csc(mm);
      for (int k=0; k<nq; k++){
	m           = k*nq;
      for (int l=0; l<nq; l++)
	bu(l,k,j,mm) = aa(l,k,j,nn);
      }
    }

    // reset working array to zero
    for (int k=0; k<nq; k++){
      m             = k*nq;
    for (int l=0; l<nq; l++)
      aa(l,k,j,c1)  = 0.;
    }
  }}

  aa.deallocate();


  // surface boundaries
  j  = 0;
  jp = 1;
  double eps=1.e-14;
  for (int n=0; n<nFaces-nGfaces; n++){
    nx[0]       =-facu(0,j,n);
    nx[1]       =-facu(1,j,n);
    ds          = sqrt(nx[0]*nx[0]+nx[1]*nx[1]);
    if (fabs(ds) < eps){// on sharp corners
      for (int k=0; k<nq; k++){
	i           = k*nq;
      for (int l=0; l<nq; l++){
	dd(l,k,j,n) = 0.;
	dp(l,k,j,n) = 0.;
      }}
      for (int k=0; k<nq; k++) dd(k,k,j,n) = 1.;
    }
    else{
      ds        = 1./ds;
      nx[0]    *= ds;
      nx[1]    *= ds;
      sys->lhsBCVectorSelfJacobian(npts,&fTag(n),&nx[0],&q(0,jp,n),
				   &qa(0,jp,n),&q(0,j,n),&qa(0,j,n),&a[0]);
      sys->lhsBCVectorInteriorJacobian(npts,&fTag(n),&nx[0],&q(0,jp,n),
				       &qa(0,jp,n),&q(0,j,n),&qa(0,j,n),&ai[0]);
      for (int k=0; k<nq; k++){
	i           = k*nq;
      for (int l=0; l<nq; l++){
	dd(l,k,j,n) = a[i+l];
	dp(l,k,j,n) = ai[i+l];
      }}
    }}


  // tip boundaries
  if (standAlone == 1){
    int eTag     = nBpatches-1; //tip tag
  for (int n=0; n<nFaces-nGfaces; n++){
    j            = fClip(n);
    jp           = j+1;
    nx[0]        = facu(0,j,n);
    nx[1]        = facu(1,j,n);
    ds           = 1./sqrt(nx[0]*nx[0]+nx[1]*nx[1]);
    nx[0]       *= ds;
    nx[1]       *= ds;
    sys->lhsBCVectorSelfJacobian(npts,&eTag,&nx[0],&q(0,j,n),
				 &qa(0,j,n),&q(0,jp,n),&qa(0,jp,n),&a[0]);
    sys->lhsBCVectorInteriorJacobian(npts,&eTag,&nx[0],&q(0,j,n),
				 &qa(0,j,n),&q(0,jp,n),&qa(0,jp,n),&ai[0]);
  for (int k=0; k<nq; k++){
    i            = k*nq;
  for (int l=0; l<nq; l++){
    dd(l,k,jp,n) = a[i+l];
    dm(l,k,jp,n) = ai[i+l];
  }}}}
}
