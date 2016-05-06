#include "StrandBlockSolver.h"


void StrandBlockSolver::specRadi()
{
  if (inviscid != 0){

    // initialize radi to zero
    int c1,c2,m,mp,jp,npts=1;
    double sr;
    for (int n=0; n<nFaces+nBedges; n++)
      for (int j=0; j<nPstr+2; j++) radi(j,n) = 0.;


    // unstructured faces
    for (int n=0; n<nEdges; n++){
      c1          = edge(0,n);
      c2          = edge(1,n);
    for (int j=1; j<nPstr+1; j++){
      sys->stepInvEigenvalue(npts,&facs(0,j,n),&xvs(j,n),&q(0,j,c1),
			     &q(0,j,c2),&qa(0,j,c1),&qa(0,j,c2),&sr);
      radi(j,c1) += sr;
      radi(j,c2) += sr;
    }}


    // structured faces
    for (int n=0; n<nFaces-nGfaces; n++){
    for (int j=0; j<fClip(n)+1; j++){
      jp          = j+1;
      sys->stepInvEigenvalue(npts,&facu(0,j,n),&xvu(j,n),&q(0,j,n),
			     &q(0,jp,n),&qa(0,j,n),&qa(0,jp,n),&sr);
      radi(j ,n) += sr;
      radi(jp,n) += sr;
    }}


    //set values at ghost cells
    for (int n=0; n<nFaces-nGfaces; n++){
      m          = fClip(n);
      mp         = m+1;
      radi(0 ,n) = radi(1,n);
      radi(mp,n) = radi(m,n);
    }
    for (int n=nEdges-nBedges; n<nEdges; n++){
      c1         = edge(0,n);
      c2         = edge(1,n);
      for (int j=1; j<nPstr+1; j++) radi(j,c2) = radi(j,c1);
    }}
}
