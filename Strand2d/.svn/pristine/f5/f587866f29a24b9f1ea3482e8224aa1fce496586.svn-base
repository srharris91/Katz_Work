#include "StrandBlockSolver.h"


void StrandBlockSolver::rhsSourceMMS()
{
  for (int n=0; n<nFaces-nGfaces; n++)
    for (int j=1; j<fClip(n)+1; j++)
      for (int k=0; k<nq; k++) r(k,j,n) -=(v(j,n)*s(k,j,n));

  int j=0;
  for (int n=0; n<nFaces-nGfaces; n++)
    for (int k=0; k<nq; k++) r(k,j,n) -= s(k,j,n);

  if (standAlone == 1){
    int j=nPstr+1;
    for (int n=0; n<nFaces-nGfaces; n++)
      for (int k=0; k<nq; k++) r(k,j,n) -= s(k,j,n);
  }

  for (int n=nFaces; n<nFaces+nBedges; n++)
    for (int j=1; j<fClip(n)+1; j++)
      for (int k=0; k<nq; k++) r(k,j,n) -= s(k,j,n);
}
