#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::prolong()
{
  // prolong corrections
  int eC,i,nC;
  double a;
  for (int n=0; n<nNode; n++){
    eC = nce(n,0);
    i  = nce(n,1);
    for(int jC=0; jC<nneC; jC++){
      nC = elemC[eC*nneC+jC];
      for (int k=0; k<nq; k++) q(n,k) += lqCF(i,jC)*(qC[nC*nq+k]-q0C[nC*nq+k]);
    }
  }


  // set additional variables
  sys->stepQAdd(nNode,
		&q(0,0),
		&qa(0,0));
}
