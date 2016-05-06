#include "StrandBlockSolver.h"


// [rhsDissipationInviscidCoarse]
void StrandBlockSolver::rhsDissipationInviscidCoarse()
{
  int c1,c2,fc,jp,npts=1;
  double ql[nq],qr[nq],fi[nq],fd[nq];


  // unstructured faces
  for (int n=0; n<nEdges; n++){
    c1             = edge(0,n);
    c2             = edge(1,n);
    fc             = fClip(c1);
    if (fClip(c2) > fc) fc = fClip(c2);
    for (int j=1; j<fc+1; j++){
      for (int k=0; k<nq; k++){
	ql[k]      = q(k,j,c1);
	qr[k]      = q(k,j,c2);
      }
      sys->rhsInvFlux(npts,&facs(0,j,n),&xvs(j,n),&ql[0],&qr[0],&fi[0]);
      sys->rhsDisFluxCoarse(npts,&facs(0,j,n),&xvs(j,n),&ql[0],&qr[0],&fd[0]);
      for (int k=0; k<nq; k++){
	r(k,j,c1) += (fi[k]-fd[k]);
	r(k,j,c2) -= (fi[k]-fd[k]);
      }}
  }


  // structured faces
  for (int n=0; n<nFaces-nGfaces; n++){
    for (int j=0; j<fClip(n)+1; j++){
      jp           = j+1;
      for (int k=0; k<nq; k++){
	ql[k]      = q(k,j ,n);
	qr[k]      = q(k,jp,n);
      }
      sys->rhsInvFlux(npts,&facu(0,j,n),&xvu(j,n),&ql[0],&qr[0],&fi[0]);
      sys->rhsDisFluxCoarse(npts,&facu(0,j,n),&xvu(j,n),&ql[0],&qr[0],&fd[0]);
      for (int k=0; k<nq; k++){
	r(k,j ,n) += (fi[k]-fd[k]);
	r(k,jp,n) -= (fi[k]-fd[k]);
      }}
  }
}
