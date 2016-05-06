/**
 * \brief
 * Implementation file for Class Tri2dFCManager. Individual methods
 * are contained in files with the name of the method itself, except for
 * static variables, constructors, destructors, and get and set methods, etc.
 *
 * \author
 * Aaron Katz
 *
 * \version
 * 1.0
 *
 * \date
 * 03-15-2013
 */


#include "Tri2dFCManager.h"


Tri2dFCManager* Tri2dFCManager::tri2dfc_manager_instance = NULL;
int Tri2dFCManager::iplotmesh = 0;
int Tri2dFCManager::nTri = 0;
int Tri2dFCManager::nNode = 0;
int Tri2dFCManager::nCompBd = 0;
int Tri2dFCManager::nEdgeBd = 0;
int Tri2dFCManager::nBlockSolvers = 0;
int Tri2dFCManager::nLevels = 0;
int Tri2dFCManager::iConvFile = 0;
int Tri2dFCManager::restartStep = 0;
int Tri2dFCManager::nOutput = 0;
int Tri2dFCManager::nSteps = 0;
int Tri2dFCManager::nPseudoSteps0 = 0;
int Tri2dFCManager::nPseudoSteps = 0;
int Tri2dFCManager::nRKStages = 0;
int Tri2dFCManager::mgCycle = 0;
int Tri2dFCManager::nq = 0;
int Tri2dFCManager::nDofs = 0;
double Tri2dFCManager::convLimit = 0.;
clock_t Tri2dFCManager::time0 = 0.;
clock_t Tri2dFCManager::timeS0 = 0.;
clock_t Tri2dFCManager::time = 0.;
Tri2dFCBlockSolver* Tri2dFCManager::t2dfcbs = NULL;
int** Tri2dFCManager::tri = NULL;


// [Tri2dFCManager]
Tri2dFCManager::Tri2dFCManager()
{
}
// [Tri2dFCManager]


// [~Tri2dFCManager]
Tri2dFCManager::~Tri2dFCManager()
{
}
// [~Tri2dFCManager]


// [createManager]
void Tri2dFCManager::createManager()
{
  if (!tri2dfc_manager_instance)
    tri2dfc_manager_instance = new Tri2dFCManager();
}
// [createManager]


// [getManager]
Tri2dFCManager* Tri2dFCManager::getManager()
{
  if (!tri2dfc_manager_instance) createManager();
  return(tri2dfc_manager_instance);
}
// [getManager]


// [freeManager]
void Tri2dFCManager::freeManager()
{
  if (tri2dfc_manager_instance) delete tri2dfc_manager_instance;
  tri2dfc_manager_instance = ((Tri2dFCManager*) NULL);
}
// [freeManager]


// get methods
const int& Tri2dFCManager::getRestartStep(){return(restartStep);}
const int& Tri2dFCManager::getNOutput(){return(nOutput);}
const int& Tri2dFCManager::getNSteps(){return(nSteps);}
const int& Tri2dFCManager::getNPseudoSteps(){return(nPseudoSteps);}
const int& Tri2dFCManager::getNPseudoSteps0(){return(nPseudoSteps0);}
