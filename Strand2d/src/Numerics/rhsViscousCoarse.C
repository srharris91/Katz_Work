#include "StrandBlockSolver.h"


void StrandBlockSolver::rhsViscousCoarse()
{
  int c1,c2,fc,jp,npts=1;
  double fv[nq];


  // unstructured faces
  for (int n=0; n<nEdges; n++){
    c1            = edge(0,n);
    c2            = edge(1,n);
    fc            = fClip(c1);
    if (fClip(c2) > fc) fc = fClip(c2);
    for (int j=1; j<fc+1; j++){
      sys->rhsVisFluxCoarse(npts,&facs(0,j,n),&v(j,c1),&v(j,c2),&q(0,j,c1),
			    &q(0,j,c2),&qa(0,j,c1),&qa(0,j,c2),&fv[0]);
      for (int k=0; k<nq; k++){
	r(k,j,c1) -= fv[k];
	r(k,j,c2) += fv[k];
      }}}


  // structured faces
  for (int n=0; n<nFaces-nGfaces; n++){
    for (int j=0; j<fClip(n)+1; j++){
      jp           = j+1;
      sys->rhsVisFluxCoarse(npts,&facu(0,j,n),&v(j,n),&v(jp,n),&q(0,j,n),
			    &q(0,jp,n),&qa(0,j,n),&qa(0,jp,n),&fv[0]);
      for (int k=0; k<nq; k++){
	r(k,j ,n) -= fv[k];
	r(k,jp,n) += fv[k];
      }}}
}
