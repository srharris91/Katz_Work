#include "Strand2dFCManager.h"


void Strand2dFCManager::takePseudoStep(const int& step,
				       const int& pseudoStep,
				       bool& converged)
{
  // shift time levels/initialize the solution at n+1 for unsteady case
  if (step > 0 && pseudoStep == 0){
    cout << "\nSolving physical time step " << step << endl;
    for (int level=0; level<nLevels; level++)
      blockSolver(0,level).shiftTime(step);
  }

  if (pseudoStep == 0) timeS0 = clock();


  // declare and initialize working variables
  int stage,pos=0,level=mgLevel(0),mode=mgMode(0),
    j,strand0,strand1,strandi,sweep,nSweep=2,nStrandNode;
  converged = false;


  // start MG cycle
  do{ 
    // prolong corrections
    if (mode >= 4) blockSolver(0,level).prolong();


    // perform Runge-Kutta iterations
    if (mode <= 4)
      for (sweep=0; sweep<nSweep; sweep++){
	nStrandNode = blockMesh(0,level).getNStrandNode();
	if (sweep % 2 == 0){
	  strand0 = 0;
	  strand1 = nStrandNode;
	  strandi = 1;
	}
	else{
	  strand0 = nStrandNode-1;
	  strand1 =-1;
	  strandi =-1;
	}
	j = strand0;
	while (j != strand1){
	  for (stage=0; stage<nRKStages; stage++)
	    blockSolver(0,level).solve(step,
				       pseudoStep,
				       sweep,
				       stage,
				       mode,
				       j);
	  j += strandi;
	}}


    // recompute residual for MG restriction, then restrict
    if ((mode == 1 || mode == 2 || mode == 4) && nLevels > 1){
      blockSolver(0,level).computeRHS(step,
				      pseudoStep);
      blockSolver(0,level+1).restrict();
    }

    // move to next level
    pos++; level = mgLevel(pos); mode = mgMode(pos);

    // check convergence if done with the entire MG cycle
    if (mode == 0) convHistory(step,
			       pseudoStep,
			       converged);
  } while(mode != 0);
}
