#include "StrandBlockSolver.h"


void StrandBlockSolver::rhsBoundaryFine()
{
  int c1,c2,m,j,jp,npts=1,eTag=nBpatches-1;//tip tag
  double wx[ndim],nx[ndim],ds,dx,dy,qe[nq],qae[nqa],eps=1.e-14;


  // unstructured boundaries
  for (int n=nEdges-nBedges; n<nEdges; n++){
    m           = n-(nEdges-nBedges);
    c1          = edge(0,n);
    c2          = edge(1,n);
    for (int j=1; j<fClip(c1)+1; j++){
      wx[0]     = nvs(0,j,m);
      wx[1]     = nvs(1,j,m);
      nx[0]     = facs(0,j,n);
      nx[1]     = facs(1,j,n);
      ds        = 1./sqrt(nx[0]*nx[0]+nx[1]*nx[1]);
      nx[0]    *= ds;
      nx[1]    *= ds;
      dx        = xc(0,j,c2)-xc(0,j,c1);
      dy        = xc(1,j,c2)-xc(1,j,c1);
      for (int k=0; k<nq; k++)
	qe[k] = q(k,j,c1)+brelax*lims(k,j,n)*(dx*qx(k,0,j,c1)+dy*qx(k,1,j,c1));
      sys->stepQAdd(npts,&qe[0],&qae[0]);
      sys->rhsBCVector(npts,&bTag(m),&nx[0],&wx[0],&qe[0],&qae[0],
		       &q(0,j,c2),&qa(0,j,c2),&r(0,j,c2));
    }}


  // surface and tip boundaries
  for (int n=0; n<nFaces-nGfaces; n++){
    j         = 0;
    jp        = 1;
    wx[0]     = nvu(0,0,n);
    wx[1]     = nvu(1,0,n);
    nx[0]     =-facu(0,j,n);
    nx[1]     =-facu(1,j,n);
    ds        = sqrt(nx[0]*nx[0]+nx[1]*nx[1]);
    if (fabs(ds) < eps){ // sharp corners
      for (int k=0; k<nq; k++) r(k,j,n) = 0.;
    }
    else{
      ds      = 1./ds;
      nx[0]  *= ds;
      nx[1]  *= ds;
      dx      = xc(0,j,n)-xc(0,jp,n);
      dy      = xc(1,j,n)-xc(1,jp,n);
      for (int k=0; k<nq; k++)
	qe[k] = q(k,jp,n)+brelax*limu(k,j,n)*(dx*qx(k,0,jp,n)+dy*qx(k,1,jp,n));
      sys->stepQAdd(npts,&qe[0],&qae[0]);
      sys->rhsBCVector(npts,&fTag(n),&nx[0],&wx[0],&qe[0],&qae[0],
		       &q(0,j,n),&qa(0,j,n),&r(0,j,n));
    }}

  if (standAlone == 1){
    for (int n=0; n<nFaces-nGfaces; n++){
      j         = fClip(n);
      jp        = j+1;
      wx[0]     = nvu(0,1,n);
      wx[1]     = nvu(1,1,n);
      nx[0]     = facu(0,j,n);
      nx[1]     = facu(1,j,n);
      ds        = 1./sqrt(nx[0]*nx[0]+nx[1]*nx[1]);
      nx[0]    *= ds;
      nx[1]    *= ds;
      dx        = xc(0,jp,n)-xc(0,j,n);
      dy        = xc(1,jp,n)-xc(1,j,n);
      for (int k=0; k<nq; k++)
	qe[k]   = q(k,j,n)+brelax*limu(k,j,n)*(dx*qx(k,0,j,n)+dy*qx(k,1,j,n));
      sys->stepQAdd(npts,&qe[0],&qae[0]);
      sys->rhsBCVector(npts,&eTag,&nx[0],&wx[0],&qe[0],&qae[0],
		       &q(0,jp,n),&qa(0,jp,n),&r(0,jp,n));
    }}
}
