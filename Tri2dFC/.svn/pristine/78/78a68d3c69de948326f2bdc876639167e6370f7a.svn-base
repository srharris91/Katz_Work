#include "Tri2dFCBlockSolver.h"
#include "lagrangePoly.h"


void Tri2dFCBlockSolver::metric()
{
  // node volumes
  v.allocate(nNode);
  volume();


  // face areas
  area.allocate(nEdge,2);
  areaBd.allocate(nEdgeBd,2);
  faceArea();


  // face areas for viscous terms
  if (order == 3){
    areaE.allocate(nElem,nee,2);
    faceAreaVis();
  }


  // set up FEM gradient coefficients
  if (order == 3){
    dxg.allocate(nElem,nne,nne,2);
    gx.allocate(npsp1,2);
    gxx.allocate(npsp1,3);
    gradSetup();
  }


  // boundary normals
  ln.allocate(nNodeBd,2);
  bNormal();
}
