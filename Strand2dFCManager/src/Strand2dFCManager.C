/**
 * \brief
 * Implementation file for Class Strand2dFCManager. Individual methods
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


#include "Strand2dFCManager.h"


Strand2dFCManager* Strand2dFCManager::Strand2dfc_manager_instance = NULL;
int Strand2dFCManager::nLevels = 0;
int Strand2dFCManager::iplotmesh = 0;
int Strand2dFCManager::nBlockSolvers = 0;
int Strand2dFCManager::nBlockMeshes = 0;
int Strand2dFCManager::iConvFile = 0;
int Strand2dFCManager::restartStep = 0;
int Strand2dFCManager::nOutput = 0;
int Strand2dFCManager::nSteps = 0;
int Strand2dFCManager::nPseudoSteps0 = 0;
int Strand2dFCManager::nPseudoSteps = 0;
int Strand2dFCManager::nRKStages = 0;
int Strand2dFCManager::mgCycle = 0;
int Strand2dFCManager::nq = 0;
int Strand2dFCManager::nDofs = 0;
double Strand2dFCManager::convLimit = 0.;
clock_t Strand2dFCManager::time0 = 0.;
clock_t Strand2dFCManager::timeS0 = 0.;
clock_t Strand2dFCManager::time = 0.;


// [Strand2dFCManager]
Strand2dFCManager::Strand2dFCManager()
{
}
// [Strand2dFCManager]


// [~Strand2dFCManager]
Strand2dFCManager::~Strand2dFCManager()
{
}
// [~Strand2dFCManager]


// [createManager]
void Strand2dFCManager::createManager()
{
  if (!Strand2dfc_manager_instance)
    Strand2dfc_manager_instance = new Strand2dFCManager();
}
// [createManager]


// [getManager]
Strand2dFCManager* Strand2dFCManager::getManager()
{
  if (!Strand2dfc_manager_instance) createManager();
  return(Strand2dfc_manager_instance);
}
// [getManager]


// [freeManager]
void Strand2dFCManager::freeManager()
{
  if (Strand2dfc_manager_instance) delete Strand2dfc_manager_instance;
  Strand2dfc_manager_instance = ((Strand2dFCManager*) NULL);
}
// [freeManager]


// [plot]
void Strand2dFCManager::plot()
{
  if (iplotmesh != 0)
    for (int level=0; level<nLevels; level++) blockMesh(0,level).plot();
}
// [plot]


// get methods
const int& Strand2dFCManager::getRestartStep(){return(restartStep);}
const int& Strand2dFCManager::getNOutput(){return(nOutput);}
const int& Strand2dFCManager::getNSteps(){return(nSteps);}
const int& Strand2dFCManager::getNPseudoSteps(){return(nPseudoSteps);}
const int& Strand2dFCManager::getNPseudoSteps0(){return(nPseudoSteps0);}
