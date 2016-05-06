#include "Tri2dFCManager.h"


void Tri2dFCManager::takePseudoStep(const int& step,
				    const int& pseudoStep,
				    bool& converged)
{
  // shift time levels/initialize the solution at n+1 for unsteady case
  if (step > 0 && pseudoStep == 0){
    cout << "\nSolving physical time step " << step << endl;
    for (int level=0; level<nLevels; level++)
      t2dfcbs[level].shiftTime(step);
  }

  if (pseudoStep == 0) timeS0 = clock();


  // declare and initialize working variables
  int stage,pos=0,level=mgLevel(0),mode=mgMode(0);
  converged = false;


  // start MG cycle
  do{ 
    // prolong corrections
    if (mode >= 4) prolong(level);


    // perform Runge-Kutta iterations
    if (mode <= 4)
      for (stage=0; stage<nRKStages; stage++)
	t2dfcbs[level].solve(step,
			     pseudoStep,
			     stage,
			     mode);


    // recompute residual for MG restriction, then restrict
    if ((mode == 1 || mode == 2 || mode == 4) && nLevels > 1){
      mode =-mode; stage = 0;
      t2dfcbs[level].computeRHS(step,
				pseudoStep,
				stage,
				mode);
      mode =-mode;
      restrict(level);
    }

    // move to next level
    pos++; level = mgLevel(pos); mode = mgMode(pos);

    // check convergence if done with the entire MG cycle
    if (mode == 0) convHistory(step,
			       pseudoStep,
			       converged);
  } while(mode != 0);
}
