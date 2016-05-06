#include "StrandBlockSolver.h"


void StrandBlockSolver::prolong(const int& mglevel,
				const int& mode)
{
  // compute coarse level corrections and inject them into the fine level
  int j=nPstrP+2,n=nFacesP+nBedgesP,npts;
  npts = j*n;
  Array3D<double> a(nq,j,n);
  for (int n=0; n<nFacesP+nBedgesP; n++)
    for (int j=0; j<nPstrP+2; j++)
      for (int k=0; k<nq; k++) a(k,j,n) = 0.;

  int m=nPstr+1,mP=nPstrP+1,i,l;
  for (int n=0; n<nFacesP+nBedgesP; n++){
    i = f2cc(n);
    if (i >= 0){
      for (int k=0; k<nq; k++){
	a(k,0 ,n) = q(k,0,i)-q0(k,0,i);
	a(k,mP,n) = q(k,m,i)-q0(k,m,i);
      }
      for (int j=1; j<(*fClipP)(n)+1; j++){
	l = f2cs(j);
	for (int k=0; k<nq; k++)
	  a(k,j,n) = q(k,l,i)-q0(k,l,i);
      }}}

  // set corner corrections to zero
  for (int n=nFacesP-nGfacesP; n<nFacesP+nBedgesP; n++)
    for (int k=0; k<nq; k++){
      a(k,0 ,n) = 0.;
      a(k,mP,n) = 0.;
    }

  // set ghost corrections to zero
  for (int n=nFacesP-nGfacesP; n<nFacesP; n++)
    for (int j=0; j<nPstrP+2; j++)
      for (int k=0; k<nq; k++)
	a(k,j,n) = 0.;




  // smooth corrections
  int smooth=0;
  if (smooth == 1){
    double fsmoothu=.25,fsmooths=.125;
    // unstructured direction
    int j=nPstrP+2,n=nFacesP+nBedgesP,c1,c2,fc,qk,gn,gr;
    Array3D<double> dqn(nq,j,n);
    Array3D<double> dqi(nq,j,n);
    Array2D<double> g  (   j,n);

    for (int n=0; n<nFacesP+nBedgesP; n++)
      for (int j=0; j<nPstrP+2; j++){
	g(j,n) = 0.;
	for (int k=0; k<nq; k++) dqn(k,j,n) = a(k,j,n);
      }

    for (int n=0; n<nEdgesP-nBedgesP; n++){
      c1         = (*edgeP)(0,n);
      c2         = (*edgeP)(1,n);
      fc         = (*fClipP)(c1);
      if ((*fClipP)(c2) > fc) fc = (*fClipP)(c2);
      for (int j=1; j<fc+1; j++){
	g(j,c1) += fsmoothu;
	g(j,c2) += fsmoothu;
      }}


    for (int l=0; l<2; l++){
      for (int n=0; n<nFacesP+nBedgesP; n++)
	for (int j=0; j<nPstrP+2; j++)
	  for (int k=0; k<nq; k++) dqi(k,j,n) = 0.;

      for (int n=0; n<nEdgesP-nBedgesP; n++){
	c1      = (*edgeP)(0,n);
	c2      = (*edgeP)(1,n);
	fc      = (*fClipP)(c1);
	if ((*fClipP)(c2) > fc) fc = (*fClipP)(c2);
	for (int j=1; j<fc+1; j++)
	  for (int k=0; k<nq; k++){
	    qk           = a(k,j,c2)-a(k,j,c1);
	    dqi(k,j,c1) += qk;
	    dqi(k,j,c2) -= qk;
	  }}

      for (int n=0; n<nFacesP-nGfacesP; n++){
	for (int j=1; j<(*fClipP)(n)+1; j++){
	  gn         = g(j,n);
	  gr         = 1./(1.+gn);
	  for (int k=0; k<nq; k++){
	    a(k,j,n) =(dqn(k,j,n)+gn*a(k,j,n)+fsmoothu*dqi(k,j,n))*gr;
	  }}}}

    dqn.deallocate();
    dqi.deallocate();
    g.deallocate();
    // strand direction
    int jm,jp;
    double eb[nPstrP+2],em,ep,ee;
    for (int n=0; n<nFacesP-nGfacesP; n++){
      em        = fsmooths;
      ep        = fsmooths;
      ee        = 1./(1.+ep);
      eb[0]     = ee*ep;
      for (int k=0; k<nq; k++) a(k,0,n) *= ee;
      for (int j=1; j<(*fClipP)(n)+2; j++){
	jm      = j-1;
	em      = ep;
	ep      = fsmooths;
	if (j == (*fClipP)(n)+1) ep = 0.;
	ee      = 1./(1.+em+ep-em*eb[jm]);
	eb[j]   = ee*ep;
	for (int k=0; k<nq; k++) a(k,j,n) = ee*(a(k,j,n)+em*a(k,jm,n));
      }
      for (int j=(*fClipP)(n); j>=0; j--){
	jp = j+1;
	for (int k=0; k<nq; k++) a(k,j,n) += (eb[j]*a(k,jp,n));
      }}}




  // apply correction
  for (int n=0; n<nFacesP+nBedgesP; n++)
    for (int j=0; j<(*fClipP)(n)+1; j++)
      //for (int n=0; n<nFacesP-nGfacesP; n++)
      //for (int j=1; j<nPstrP+1; j++)
      for (int k=0; k<nq; k++)
	(*qP)(k,j,n) += a(k,j,n);

  sys->stepQAdd(npts,
		&((*qP)(0,0,0)),
		&((*qaP)(0,0,0)));

  // if correcting finest level, reset compute flags
  if (mglevel == 1){
    *nodalQFlagP  = 0;
    *nodalQaFlagP = 0;
    *gradQFlagP   = 0;
    *gradQaFlagP  = 0;
    *limFlagP     = 0;
  }

  a.deallocate();
}
