#include "Tri2dFCManager.h"


void Tri2dFCManager::output(const int& step)
{
  t2dfcbs[0].output(nBlockSolvers,
		    step);
}
