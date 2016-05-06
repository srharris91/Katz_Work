#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::restrict()
{
  // restrict solution
  int eF,i,nF;
  q.set(0.);
  for (int n=0; n<nNode; n++){
    eF = nfe(n,0);
    i  = nfe(n,1);
    for(int jF=0; jF<nneF; jF++){
      nF = elemF[eF*nneF+jF];
      for (int k=0; k<nq; k++) q(n,k) += lqFC(i,jF)*qF[nF*nq+k];
    }}


  // pick off solution at first descent into a coarse level
  for(int n=0; n<nNode; n++)
    for (int k=0; k<nq; k++) q0(n,k) = q(n,k);


  // set additional variables
  sys->stepQAdd(nNode,
		&q(0,0),
		&qa(0,0));


  // restrict residual
  fwc.set(0.);
  double a;
  for (int n=0; n<nNodeF; n++)
    for (int j=0; j<nfn[n]; j++){
      i = nfn1[n][j];
      a = nfn2[n][j];
      for (int k=0; k<nq; k++) fwc(i,k) += a*rF[n*nq+k];
    }
}
