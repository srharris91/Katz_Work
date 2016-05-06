#include "StrandBlockSolver.h"


void StrandBlockSolver::nodalQa(const int& mglevel)
{
  if (nodalQaFlag == 0 && mglevel == 0){
    int nn,mm,jj,il;
    double qamin[nqa],qamax[nqa],a=1.,b=1.;//a=gradClip,b=1./gradClip;
    for (int n=0; n<nNodes-nGnodes; n++){
      mm = ncsp(n);
    for (int j=0; j<nPstr+1; j++){
      if      (j == 0    ) jj = 1;
      else if (j == nPstr) jj = nPstr-1;
      else                 jj = j;
    for (int k=0; k<nqa; k++) qap(k,j,n) = 0.;
    for (int k=0; k<nqa; k++) qamax[k]   =-1.e14;
    for (int k=0; k<nqa; k++) qamin[k]   = 1.e14;
    for (int i=0; i<2; i++){
      indlsp(i,j,n,il);
    for (int m=0; m<mm; m++){
      nn = csp[n][m];
      for (int k=0; k<nqa; k++) qap(k,j,n) += lsp[il][m]*qa(k,jj,nn);
      for (int k=0; k<nqa; k++) if (qa(k,jj,nn) > qamax[k]) qamax[k] = qa(k,jj,nn);
      for (int k=0; k<nqa; k++) if (qa(k,jj,nn) < qamin[k]) qamin[k] = qa(k,jj,nn);
    }
    jj++;
    }
    for (int k=0; k<nqa; k++) if (qap(k,j,n) > a*qamax[k]) qap(k,j,n) = qamax[k];
    for (int k=0; k<nqa; k++) if (qap(k,j,n) < b*qamin[k]) qap(k,j,n) = qamin[k];
    }}


    // average sharp corner values
    if (nSharp > 0){
      Array2D<double> qapv(nqa,nSharp);
      Array1D<double> deg(   nSharp);
      for (int n=0; n<nSharp; n++)
	for (int k=0; k<nqa; k++) qapv(k,n) = 0.;
      for (int n=0; n<nSharp; n++) deg(n) = 0.;
      int j=0;
      int i,im;
      for (int n=0; n<nNodes; n++){
	i = sFlag(n);
	if (i > 0){
	  im = i-1;
	  for (int k=0; k<nqa; k++) qapv(k,im) += qap(k,j,n);
	  deg(im) += 1.;
	}}
      for (int n=0; n<nNodes; n++){
	i = sFlag(n);
	if (i > 0){
	  im = i-1;
	  for (int k=0; k<nqa; k++) qap(k,j,n) = qapv(k,im)/deg(im);
	}}
      qapv.deallocate();
      deg.deallocate();
    }

    int n1,n2;
    for (int n=0; n<nFaces-nGfaces; n++){
      n1 = face(0,n);
      n2 = face(1,n);
      if (sFlag(n1) > 0 && sFlag(n2) > 0)
	for (int k=0; k<nqa; k++) qa(k,0,n) = qap(k,0,n1);
    }
  }
  nodalQaFlag = 1;



  /*
  if (nodalQaFlag == 0 && mglevel == 0){
    int nn,mm,jj,il;
    for (int n=0; n<nNodes-nGnodes; n++){
      mm = ncsp(n);
    for (int j=0; j<nPstr+1; j++){
      if      (j == 0    ) jj = 1;
      else if (j == nPstr) jj = nPstr-1;
      else                 jj = j;
    for (int k=0; k<nqa; k++) qap(k,j,n) = 0.;
    for (int i=0; i<2; i++){
      indlsp(i,j,n,il);
    for (int m=0; m<mm; m++){
      nn = csp[n][m];
      for (int k=0; k<nqa; k++) qap(k,j,n) += lsp[il][m]*qa(k,jj,nn);
    }
    jj++;
    }}}

    // average sharp corner values
    if (nSharp > 0){
      Array2D<double> qapv(nqa,nSharp);
      Array1D<double> deg(   nSharp);
      for (int n=0; n<nSharp; n++)
	for (int k=0; k<nqa; k++) qapv(k,n) = 0.;
      for (int n=0; n<nSharp; n++) deg(n) = 0.;
      int j=0;
      int i,im;
      for (int n=0; n<nNodes; n++){
	i = sFlag(n);
	if (i > 0){
	  im = i-1;
	  for (int k=0; k<nqa; k++) qapv(k,im) += qap(k,j,n);
	  deg(im) += 1.;
	}}
      for (int n=0; n<nNodes; n++){
	i = sFlag(n);
	if (i > 0){
	  im = i-1;
	  for (int k=0; k<nqa; k++) qap(k,j,n) = qapv(k,im)/deg(im);
	}}
      qapv.deallocate();
      deg.deallocate();
    }
  }
  nodalQaFlag = 1;
  */
}
