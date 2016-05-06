#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::connectivity()
{
  // extract edges from the surface elements
  edgeExtract();


  // put boundary nodes at the end, and update all data structures
  reorderNodes();


  // form inter-element gradient stencil
  if (surfOrder > 1) gradStencil();

  // form source term stencil
  sourceStencil();
}
