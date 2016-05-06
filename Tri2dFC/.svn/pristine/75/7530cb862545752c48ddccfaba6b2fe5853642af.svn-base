#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::solve(const int& step,
			       const int& pseudoStep,
			       const int& stage,
			       const int& mode)
{
  // compute RHS
  computeRHS(step,
	     pseudoStep,
	     stage,
	     mode);


  // save first stage q, on the fine level save q for RMS purposes
  if (stage == 0){
    for(int n=0; n<nNode; n++)
      for (int k=0; k<nq; k++) qn(n,k) = q(n,k);
    //if (level == 0)
      for(int n=0; n<nNode; n++)
	for (int k=0; k<nq; k++) q0(n,k) = q(n,k);
  }


  // compute local pseudo-time step divided by volume using spectral radius
  if (stage == 0) pseudoTime();


  // multiply the RHS by the time step and update
  update(step,
	 stage);
}
