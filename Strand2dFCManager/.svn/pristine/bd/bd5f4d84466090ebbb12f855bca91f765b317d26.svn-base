#include "Strand2dFCManager.h"


void Strand2dFCManager::finalize()
{
  for (int level=0; level<nLevels; level++) blockMesh(0,level).finalize();
  blockMesh.deallocate();
  mgLevel.deallocate();
  mgMode.deallocate();
}
