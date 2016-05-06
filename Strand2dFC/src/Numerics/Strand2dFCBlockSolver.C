/**
 * \brief
 * Implementation file for Class Strand2dFCBlockSolver. Individual methods
 * are contained in files with the name of the method itself, except for
 * static variables, constructors, destructors, and get and set methods.
 *
 * \author
 * Aaron Katz
 *
 * \version
 * 1.0
 *
 * \date
 * 2013-01-01
 */


#include "Strand2dFCBlockSolver.h"


int    Strand2dFCBlockSolver::iPrint=0;
int    Strand2dFCBlockSolver::iTest=0;
int    Strand2dFCBlockSolver::iDebug=0;
int    Strand2dFCBlockSolver::iConvFile=0;
int    Strand2dFCBlockSolver::iSolnFile=0;
int    Strand2dFCBlockSolver::iResdFile=0;
int    Strand2dFCBlockSolver::iErrFile=0;
int    Strand2dFCBlockSolver::iSurfFile=0;
int    Strand2dFCBlockSolver::standAlone=0;
int    Strand2dFCBlockSolver::restartStep=0;
int    Strand2dFCBlockSolver::nRestart=0;
int    Strand2dFCBlockSolver::nOutput=0;
int    Strand2dFCBlockSolver::nSteps=0;
int    Strand2dFCBlockSolver::nPseudoSteps=0;
int    Strand2dFCBlockSolver::nPseudoSteps0=0;
int    Strand2dFCBlockSolver::nLinearSteps=0;
int    Strand2dFCBlockSolver::nRKStages=0;
int    Strand2dFCBlockSolver::implicit=0;
int    Strand2dFCBlockSolver::gradMethod=0;
int    Strand2dFCBlockSolver::mgCycle=0;
int    Strand2dFCBlockSolver::limiter=0;
int    Strand2dFCBlockSolver::timeAcc=0;
double Strand2dFCBlockSolver::dtUnsteady=0.;
double Strand2dFCBlockSolver::cfl=0.;
double Strand2dFCBlockSolver::vnn=0.;
double Strand2dFCBlockSolver::smooth=0.;
double Strand2dFCBlockSolver::convLimit=0.;
double Strand2dFCBlockSolver::relax=0.;


// [Strand2dFCBlockSolver]
Strand2dFCBlockSolver::Strand2dFCBlockSolver()
{
  dataInit();
}
// [Strand2dFCBlockSolver]


// [~Strand2dFCBlockSolver]
Strand2dFCBlockSolver::~Strand2dFCBlockSolver()
{
}
// [~Strand2dFCBlockSolver]


// get methods
const int& Strand2dFCBlockSolver::getIConvFile(){return(iConvFile);}
const int& Strand2dFCBlockSolver::getRestartStep(){return(restartStep);}
const int& Strand2dFCBlockSolver::getNOutput(){return(nOutput);}
const int& Strand2dFCBlockSolver::getNSteps(){return(nSteps);}
const int& Strand2dFCBlockSolver::getNPseudoSteps(){return(nPseudoSteps);}
const int& Strand2dFCBlockSolver::getNPseudoSteps0(){return(nPseudoSteps0);}
const int& Strand2dFCBlockSolver::getNRKStages(){return(nRKStages);}
const int& Strand2dFCBlockSolver::getMgCycle(){return(mgCycle);}
const int& Strand2dFCBlockSolver::getNq(){return(nq);}
const int& Strand2dFCBlockSolver::getNDofs(){return(nDofs);}
const int& Strand2dFCBlockSolver::getMeshOrder0(){return(meshOrder0);}
const double& Strand2dFCBlockSolver::getConvLimit(){return(convLimit);}
int* Strand2dFCBlockSolver::getSurfElem0(){return(&surfElem0(0,0));}
int* Strand2dFCBlockSolver::getBndNode(){return(&bndNode(0));}
int* Strand2dFCBlockSolver::getClip(){return(&clip(0));}
double* Strand2dFCBlockSolver::getRms(){
  for (int k=0; k<nq; k++) rms(k) = 0.;
  for(int n=0; n<nSurfNode; n++)
    for(int j=0; j<nStrandNode; j++)
      for (int k=0; k<nq; k++) rms(k) += pow((q(n,j,k)-q0(n,j,k)),2);
  for (int k=0; k<nq; k++) rms(k) /= pow(rmsNorm(k),2);
  return(&rms(0));
}
double* Strand2dFCBlockSolver::getQ(){return(&q(0,0,0));}
double* Strand2dFCBlockSolver::getQ0(){return(&q0(0,0,0));}
double* Strand2dFCBlockSolver::getR(){return(&r(0,0,0));}
