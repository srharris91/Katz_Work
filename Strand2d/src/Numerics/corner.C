#include "StrandBlockSolver.h"


void StrandBlockSolver::corner()
{
  // find sharp corners surrounding each cell
  int i,j=nFaces+nBedges;
  nssc0.allocate(j);
  ssc0 = new int*[j];
  for (int n=0; n<nFaces+nBedges; n++) nssc0(n) = 0.;
  int flag[nSharp];
  for (int n=0; n<nSharp; n++) flag[n] =-1;
  for (int n=0; n<nFaces+nBedges; n++){
    for (int m=0; m<npsc(n); m++){
      i = psc[n][m];
      j = sFlag(i);
      if (j > 0 && flag[j-1] != n){
	nssc0(n) += 1;
	flag[j-1] = n;
      }}}
  for (int n=0; n<nFaces+nBedges; n++){
    if (nssc0(n) > 0) ssc0[n] = new int[nssc0(n)];
    else ssc0[n] = NULL;
  }
  for (int n=0; n<nFaces+nBedges; n++) nssc0(n) = 0.;
  for (int n=0; n<nSharp; n++) flag[n] =-1;
  for (int n=0; n<nFaces+nBedges; n++){
    for (int m=0; m<npsc(n); m++){
      i = psc[n][m];
      j = sFlag(i);
      if (j > 0 && flag[j-1] != n){
	ssc0[n][nssc0(n)] = j-1;
	nssc0(n)         += 1;
	flag[j-1]         = n;
      }}}
  //for (int n=0; n<nFaces+nBedges; n++){
  //  cout << n;
  //  for (int m=0; m<nssc0(n); m++) cout << "\t" << ssc0[n][m];
  //  cout << endl;
  //}
  //exit(0);


  // find nodes surrounding each sharp corner
  npsc0.allocate(nSharp);
  psc0 = new int*[nSharp];
  for (int n=0; n<nSharp; n++) npsc0(n) = 0;
  for (int n=0; n<nNodes; n++){
    j = sFlag(n)-1;
    if (j >= 0) npsc0(j) += 1;
  }
  for (int n=0; n<nSharp; n++) psc0[n] = new int[npsc0(n)];
  for (int n=0; n<nSharp; n++) npsc0(n) = 0.;
  for (int n=0; n<nNodes; n++){
    j = sFlag(n)-1;
    if (j >= 0){
      psc0[j][npsc0(j)] = n;
      npsc0(j)         += 1;
    }}
  //for (int n=0; n<nSharp; n++){
  //  cout << n;
  //  for (int m=0; m<npsc0(n); m++) cout << "\t" << psc0[n][m];
  //  cout << endl;
  //}
  //exit(0);


  // find cells surrounding each sharp corner
  int n1,n2;
  ncsc0.allocate(nSharp);
  csc0 = new int*[nSharp];
  for (int n=0; n<nSharp; n++) ncsc0(n) = 0;
  for (int n=0; n<nSharp; n++) flag[n] =-1;
  for (int n=0; n<nFaces-nGfaces; n++){
    n1 = face(0,n);
    n2 = face(1,n);
    j  = sFlag(n1)-1;
    if (j >= 0 && flag[j] != n){
      ncsc0(j) += 1;
      flag[j]   = n;
    }
    j  = sFlag(n2)-1;
    if (j >= 0 && flag[j] != n){
      ncsc0(j) += 1;
      flag[j]   = n;
    }}
  for (int n=0; n<nSharp; n++) csc0[n] = new int[ncsc0(n)];
  for (int n=0; n<nSharp; n++) ncsc0(n) = 0.;
  for (int n=0; n<nSharp; n++) flag[n] =-1;
  for (int n=0; n<nFaces-nGfaces; n++){
    n1 = face(0,n);
    n2 = face(1,n);
    j  = sFlag(n1)-1;
    if (j >= 0 && flag[j] != n){
      csc0[j][ncsc0(j)] = n;
      ncsc0(j)         += 1;
      flag[j]          = n;
    }
    j  = sFlag(n2)-1;
    if (j >= 0 && flag[j] != n){
      csc0[j][ncsc0(j)] = n;
      ncsc0(j)         += 1;
      flag[j]           = n;
    }}
  //for (int n=0; n<nSharp; n++){
  //  cout << n;
  //  for (int m=0; m<ncsc0(n); m++) cout << "\t" << csc0[n][m];
  //  cout << endl;
  //}
}
