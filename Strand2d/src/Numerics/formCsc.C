#include "StrandBlockSolver.h"


void StrandBlockSolver::formCsc()
{
  int i,j,k,c1,c2;
  for (int n=0; n<nFaces+nBedges+1; n++) ncsc(n) = 0;
  for (int n=0; n<nEdges; n++)
    for (int k=0; k<2; k++){
      j        = edge(k,n)+1;
      ncsc(j) += 1;
    }
  for (int n=1; n<nFaces+nBedges+1; n++) ncsc(n) += ncsc(n-1);
  i = nFaces+nBedges;
  j = ncsc(i);
  k = nPstr+2;
  csc.allocate(j);
  bu.allocate(nq,nq,k,j);
  for (int n=0; n<nEdges; n++){
    c1        = edge(0,n);
    c2        = edge(1,n);
    i         = ncsc(c1);
    csc(i)    = c2;
    ncsc(c1) += 1;
    i         = ncsc(c2);
    csc(i)    = c1;
    ncsc(c2) += 1;
  }
  for (int n=nFaces+nBedges; n>0; n--) ncsc(n) = ncsc(n-1);
  ncsc(0) = 0;
}
