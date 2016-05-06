#include "StrandBlockSolver.h"


void StrandBlockSolver::nodalQ(const int& mglevel)
{
  if (nodalQFlag == 0 && mglevel == 0){
    int nn,mm,jj,il;
    double qmin[nq],qmax[nq],a=1.,b=1.;//a=gradClip,b=1./gradClip;
    for (int n=0; n<nNodes-nGnodes; n++){
      mm = ncsp(n);
    for (int j=0; j<nPstr+1; j++){
      if      (j == 0    ) jj = 1;
      else if (j == nPstr) jj = nPstr-1;
      else                 jj = j;
    for (int k=0; k<nq; k++) qp(k,j,n) = 0.;
    for (int k=0; k<nq; k++) qmax[k] =-1.e14;
    for (int k=0; k<nq; k++) qmin[k] = 1.e14;
    for (int i=0; i<2; i++){
      indlsp(i,j,n,il);
    for (int m=0; m<mm; m++){
      nn = csp[n][m];
      for (int k=0; k<nq; k++) qp(k,j,n) += lsp[il][m]*q(k,jj,nn);
      for (int k=0; k<nq; k++) if (q(k,jj,nn) > qmax[k]) qmax[k] = q(k,jj,nn);
      for (int k=0; k<nq; k++) if (q(k,jj,nn) < qmin[k]) qmin[k] = q(k,jj,nn);
    }
    jj++;
    }
    for (int k=0; k<nq; k++) if (qp(k,j,n) > a*qmax[k]) qp(k,j,n) = qmax[k];
    for (int k=0; k<nq; k++) if (qp(k,j,n) < b*qmin[k]) qp(k,j,n) = qmin[k];
    }}


    // average sharp corner values
    if (nSharp > 0){
      Array2D<double> qpv(nq,nSharp);
      Array1D<double> deg(   nSharp);
      for (int n=0; n<nSharp; n++)
	for (int k=0; k<nq; k++) qpv(k,n) = 0.;
      for (int n=0; n<nSharp; n++) deg(n) = 0.;
      int j=0;
      int i,im;
      for (int n=0; n<nNodes; n++){
	i = sFlag(n);
	if (i > 0){
	  im = i-1;
	  for (int k=0; k<nq; k++) qpv(k,im) += qp(k,j,n);
	  deg(im) += 1.;
	}}
      for (int n=0; n<nNodes; n++){
	i = sFlag(n);
	if (i > 0){
	  im = i-1;
	  for (int k=0; k<nq; k++) qp(k,j,n) = qpv(k,im)/deg(im);
	}}
      qpv.deallocate();
      deg.deallocate();
    }

    int n1,n2;
    for (int n=0; n<nFaces-nGfaces; n++){
      n1 = face(0,n);
      n2 = face(1,n);
      if (sFlag(n1) > 0 && sFlag(n2) > 0)
	for (int k=0; k<nq; k++) q(k,0,n) = qp(k,0,n1);
    }
  }
  nodalQFlag = 1;
}
