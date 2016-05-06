#include "Tri2dFCManager.h"


void Tri2dFCManager::finalize()
{
  if (t2dfcbs){
    for (int n=0; n<nBlockSolvers*nLevels; n++) t2dfcbs[n].finalize();
    delete [] t2dfcbs;
    t2dfcbs = NULL;
  }
  for (int n=0; n<nTri; n++){
    if (tri[n]) delete [] tri[n];
    tri[n] = NULL;
  }
  delete [] tri;
  tri = NULL;
  x.deallocate();
  edgeBd.deallocate();
  orders.deallocate();
  mgLevel.deallocate();
  mgMode.deallocate();
}
