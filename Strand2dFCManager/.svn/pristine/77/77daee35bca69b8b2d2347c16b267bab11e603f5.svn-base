#include "Strand2dFCManager.h"


void Strand2dFCManager::finalizeSolver()
{
  for (int n=0; n<nBlockSolvers; n++)
    for (int l=0; l<nLevels; l++) blockSolver(n,l).finalize();
  blockSolver.deallocate();
}
