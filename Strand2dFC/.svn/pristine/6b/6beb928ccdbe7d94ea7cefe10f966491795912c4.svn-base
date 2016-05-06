#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::restrict()
{
  // I'm on the coarse level here, grabbing data from the fine level
  // restrict solution
  int nF,nm,jF,jm,nStrandNodeF=2*(nStrandNode-1)+1;
  double wF;
  q.set(0.);
  for (int n=0; n<nSurfNode; n++)
    for (int i=psp2MGS(n); i<psp2MGS(n+1); i++){
      nF = psp1MGS(i);
      wF = wsp1MGS(i);
      nm = nF*nStrandNodeF*nq;
      for (int j=0; j<nStrandNode; j++){
	jF = 2*j;
	jm = jF*nq;
	for (int k=0; k<nq; k++) q(n,j,k) += wF*qF[nm+jm+k];
      }}


  // pick off solution at first descent into a coarse level
  for(int n=0; n<nSurfNode; n++)
    for (int j=0; j<nStrandNode; j++)
      for (int k=0; k<nq; k++) q0(n,j,k) = q(n,j,k);


  // set additional variables
  sys->stepQAdd(nSurfNode*nStrandNode,
		&q(0,0,0),
		&qa(0,0,0));


  // restrict residual
  int jmm,jmp,nSurfNodeF=nSurfElem0*meshOrder0F+nBndNode/2;
  Array3D<double> a(nSurfNodeF,nStrandNode,nq);
  for (int n=0; n<nSurfNodeF; n++){
    nm = n*nStrandNodeF*nq;

    jF = 0;
    jm = jF*nq;
    for (int k=0; k<nq; k++) a(n,0,k) = 0.;//rF[nm+jm+k];

    for (int j=1; j<nStrandNode-1; j++){
      jF  = 2*j;
      jm  = jF*nq;
      jmm =(jF-1)*nq;
      jmp =(jF+1)*nq;
      for (int k=0; k<nq; k++)
	a(n,j,k) = .5*rF[nm+jm+k]+.25*(rF[nm+jmm+k]+rF[nm+jmp+k]);
    }

    jF  = 2*(nStrandNode-1);
    jm  = jF*nq;
    for (int k=0; k<nq; k++) a(n,nStrandNode-1,k) = 0.;//rF[nm+jm+k];
  }

  int nC;
  double wC;
  fwc.set(0.);
  for (int n=0; n<nSurfNodeF; n++)
    for (int i=psp2MGR(n); i<psp2MGR(n+1); i++){
      nC = psp1MGR(i);
      wC = wsp1MGR(i);
      for (int j=0; j<nStrandNode; j++)
	for (int k=0; k<nq; k++) fwc(nC,j,k) += wC*a(n,j,k);
    }


  // clean up
  a.deallocate();
}
