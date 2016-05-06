#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::prolong()
{
  // I'm on the fine level here, grabbing data from the coarse level
  // prolong corrections
  int nm,jF,jm,jp,
    nSurfNodeC=nSurfElem0*meshOrder0C+nBndNode/2,
    nStrandNodeC=(nStrandNode-1)/2+1;
  double xF,xm,xp,wm,wp,dqm,dqp;
  Array3D<double> a(nSurfNodeC,nStrandNode,nq);
  for (int n=0; n<nSurfNodeC; n++){
    nm = n*nStrandNodeC*nq;
    for (int j=0; j<nStrandNodeC; j++){ //coincident nodes
      jF = j*2;
      jm = j*nq;
      for (int k=0; k<nq; k++) a(n,jF,k) = qC[nm+jm+k]-q0C[nm+jm+k];
    }
    for (int j=0; j<nStrandNodeC-1; j++){ //mid nodes
      jF = j*2;
      jm = j*nq;
      jp =(j+1)*nq;
      xF = strandX(jF+1);
      xm = strandX(jF  );
      xp = strandX(jF+2);
      wm =(xp-xF)/(xp-xm);
      wp =(xF-xm)/(xp-xm);
      for (int k=0; k<nq; k++){
	dqm         = qC[nm+jm+k]-q0C[nm+jm+k];
	dqp         = qC[nm+jp+k]-q0C[nm+jp+k];
	a(n,jF+1,k) = wm*dqm+wp*dqp;
      }}}

  int nC;
  double wC;
  for (int n=0; n<nSurfNode; n++)
    for (int i=psp2MGC(n); i<psp2MGC(n+1); i++){
      nC = psp1MGC(i);
      wC = wsp1MGC(i);
      for (int j=0; j<nStrandNode; j++)
	for (int k=0; k<nq; k++) q(n,j,k) += wC*a(nC,j,k);
    }


  // set additional variables
  sys->stepQAdd(nSurfNode*nStrandNode,
		&q(0,0,0),
		&qa(0,0,0));


  // clean up
  a.deallocate();
}
