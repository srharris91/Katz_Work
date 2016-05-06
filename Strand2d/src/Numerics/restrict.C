#include "StrandBlockSolver.h"


void StrandBlockSolver::restrict(const int& mglevel,
				 const int& mode)
{
  // restrict solution with volume weighting
  for (int n=0; n<nFaces+nBedges; n++)
    for (int j=0; j<nPstr+2; j++)
      for (int k=0; k<nq; k++) q(k,j,n) = 0.;

  int i,l,m=nPstr+1,mP=nPstrP+1;
  for (int n=0; n<nFacesP+nBedgesP; n++){
    i = f2cc(n);
  if (i >= 0){
  for (int k=0; k<nq; k++){
    q(k,0,i) += ((*qP)(k,0 ,n)*(*vP)(0 ,n));
    q(k,m,i) += ((*qP)(k,mP,n)*(*vP)(mP,n));
  }
  for (int j=1; j<nPstrP+1; j++){
    l = f2cs(j);
  for (int k=0; k<nq; k++)
    q(k,l,i) += ((*qP)(k,j,n)*(*vP)(j,n));
  }}}

  double vr;
  for (int n=0; n<nFaces+nBedges; n++)
    for (int j=0; j<nPstr+2; j++){
      vr = 1./v(j,n);
      for (int k=0; k<nq; k++) q(k,j,n) *= vr;
    }

  int npts=(nPstr+2)*(nFaces+nBedges);
  sys->stepQAdd(npts,
		&q(0,0,0),
		&qa(0,0,0));


  // pick off solution at first descent into a coarse level
  for (int n=0; n<nFaces+nBedges; n++)
    for (int j=0; j<nPstr+2; j++)
      for (int k=0; k<nq; k++) q0(k,j,n) = q(k,j,n);


  // restrict non-linear interior residuals via agglomeration
  for (int n=0; n<nFaces+nBedges; n++)
    for (int j=0; j<nPstr+2; j++)
      for (int k=0; k<nq; k++) fwc(k,j,n) = 0.;

  for (int n=0; n<nFacesP-nGfacesP; n++){
    i = f2cc(n);
    //for (int j=1; j<nPstrP+1; j++){
    for (int j=1; j<(*fClipP)(n)+1; j++){
      l = f2cs(j);
      for (int k=0; k<nq; k++)
	fwc(k,l,i) += ((*rP)(k,j,n));
    }}


  // restrict non-linear boundary residuals via volume weighting
  for (int n=nFacesP; n<nFacesP+nBedgesP; n++){
    i = f2cc(n);
    //for (int j=1; j<nPstrP+1; j++){
    for (int j=1; j<(*fClipP)(n)+1; j++){
      l = f2cs(j);
      for (int k=0; k<nq; k++)
	fwc(k,l,i) += ((*rP)(k,j,n)*(*vP)(j,n));
    }}
  for (int n=nFaces; n<nFaces+nBedges; n++)
    //for (int j=1; j<nPstr+1; j++){
    for (int j=1; j<fClip(n)+1; j++){
      vr = 1./v(j,n);
      for (int k=0; k<nq; k++) fwc(k,j,n) *= vr;
    }

  for (int n=0; n<nFacesP-nGfacesP; n++){
    i = f2cc(n);
    for (int k=0; k<nq; k++){
      fwc(k,0,i) += ((*rP)(k,0 ,n)*(*vP)(0 ,n));
    }}
  double vr0,vrm;
  for (int n=0; n<nFaces-nGfaces; n++){
    vr0 = v(0,n);
    vrm = v(m,n);
    for (int k=0; k<nq; k++){
      fwc(k,0,n) *= vr0;
    }}

  if (standAlone == 1){
    for (int n=0; n<nFacesP-nGfacesP; n++){
      i = f2cc(n);
      for (int k=0; k<nq; k++){
	fwc(k,m,i) += ((*rP)(k,mP,n)*(*vP)(mP,n));
      }}
    double vr0,vrm;
    for (int n=0; n<nFaces-nGfaces; n++){
      vrm = v(m,n);
      for (int k=0; k<nq; k++){
	fwc(k,m,n) *= vrm;
      }}
  }

  // set corner and ghost cells to 0.
  for (int n=nFaces-nGfaces; n<nFaces+nBedges; n++)
    for (int k=0; k<nq; k++){
      fwc(k,0,n) = 0.;
      fwc(k,m,n) = 0.;
    }
}
