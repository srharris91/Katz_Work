#include "StrandBlockSolver.h"


void StrandBlockSolver::specRadv()
{
  if (viscous != 0){

    // initialize radv to zero
    int c1,c2,m,mp,jp,npts=1;
    double sr;
    for (int n=0; n<nFaces+nBedges; n++)
      for (int j=0; j<nPstr+2; j++) radv(j,n) = 0.;


    // unstructured faces
    for (int n=0; n<nEdges; n++){
      c1          = edge(0,n);
      c2          = edge(1,n);
    for (int j=1; j<nPstr+1; j++){
      sys->stepVisEigenvalue(npts,&facs(0,j,n),&q(0,j,c1),&q(0,j,c2),
			     &qa(0,j,c1),&qa(0,j,c2),&v(j,c1),&v(j,c2),&sr);
      radv(j,c1) += sr;
      radv(j,c2) += sr;
    }}


    // structured faces
    for (int n=0; n<nFaces-nGfaces; n++){
    for (int j=0; j<fClip(n)+1; j++){
      jp          = j+1;
      sys->stepVisEigenvalue(npts,&facu(0,j,n),&q(0,j,n),&q(0,jp,n),
			     &qa(0,j,n),&qa(0,jp,n),&v(j,n),&v(jp,n),&sr);
      radv(j ,n) += sr;
      radv(jp,n) += sr;
    }}


    //set values at ghost cells
    for (int n=0; n<nFaces-nGfaces; n++){
      m          = fClip(n);
      mp         = m+1;
      radv(0 ,n) = radv(1,n);
      radv(mp,n) = radv(m,n);
    }
    for (int n=nEdges-nBedges; n<nEdges; n++){
      c1         = edge(0,n);
      c2         = edge(1,n);
      for (int j=1; j<nPstr+1; j++) radv(j,c2) = radv(j,c1);
    }}
}
