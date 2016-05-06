#include "Tri2dFCManager.h"


void Tri2dFCManager::restrict(const int& level)
{
  t2dfcbs[level+1].restrict();
}
