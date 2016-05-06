#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::metric()
{
  // 1d finite difference coefficients for the strand direction
  fdCoeff();


  // mapping terms, including nodal "volumes"
  mapping();


  // boundary normals
  bNormal();


  // gradient coefficients
  if (surfOrder > 1) gradSetup();

  // source coefficients
  sourceSetup();
}
