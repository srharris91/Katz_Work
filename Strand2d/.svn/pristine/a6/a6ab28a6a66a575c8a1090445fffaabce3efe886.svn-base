/**
 * \brief
 * Implementation file for Class StrandBlockSolver. Individual methods
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
 * 2012-08-30
 */


#include "StrandBlockSolver.h"


int    StrandBlockSolver::iPrint = 0;
int    StrandBlockSolver::iTest = 0;
int    StrandBlockSolver::iDebug = 0;
int    StrandBlockSolver::iConvFile = 0;
int    StrandBlockSolver::iSolnFile = 0;
int    StrandBlockSolver::iResdFile = 0;
int    StrandBlockSolver::iErrFile = 0;
int    StrandBlockSolver::iSurfFile = 0;
int    StrandBlockSolver::standAlone = 0;
int    StrandBlockSolver::restartStep = 0;
int    StrandBlockSolver::nRestart = 0;
int    StrandBlockSolver::nOutput = 0;
int    StrandBlockSolver::nSteps = 0;
int    StrandBlockSolver::nPseudoSteps = 0;
int    StrandBlockSolver::nPseudoSteps0 = 0;
int    StrandBlockSolver::nLinearSteps = 0;
int    StrandBlockSolver::nLevels = 0;
int    StrandBlockSolver::mgCycle = 0;
int    StrandBlockSolver::gradient = 0;
int    StrandBlockSolver::nodeVal = 0;
int    StrandBlockSolver::perturb = 0;
int    StrandBlockSolver::limiter = 0;
int    StrandBlockSolver::nRamp = 0;
double StrandBlockSolver::brelax = 0.;
double StrandBlockSolver::gradClip = 0.;
double StrandBlockSolver::dtUnsteady = 0.;
double StrandBlockSolver::cflLinear = 0.;
double StrandBlockSolver::cfl0 = 0.;
double StrandBlockSolver::vnn0 = 0.;
double StrandBlockSolver::cfl = 0.;
double StrandBlockSolver::vnn = 0.;
double StrandBlockSolver::convLimit = 0.;
double StrandBlockSolver::relax = 0.;
double StrandBlockSolver::coarseDis = 0.;


// [StrandBlockSolver]
StrandBlockSolver::StrandBlockSolver()
{
  dataInit();
}
// [StrandBlockSolver]


// [StrandBlockSolver]
StrandBlockSolver::StrandBlockSolver(const StrandBlockSolver& sbs)
{
  dataInit();
}
// [StrandBlockSolver]


// [~StrandBlockSolver]
StrandBlockSolver::~StrandBlockSolver()
{
}
// [~StrandBlockSolver]


const int& StrandBlockSolver::getIConvFile()
{
  return(iConvFile);
}

const int& StrandBlockSolver::getRestartStep()
{
  return(restartStep);
}

const int& StrandBlockSolver::getNSteps()
{
  return(nSteps);
}

const int& StrandBlockSolver::getNPseudoSteps()
{
  return(nPseudoSteps);
}

const int& StrandBlockSolver::getNPseudoSteps0()
{
  return(nPseudoSteps0);
}

const int& StrandBlockSolver::getNLinearSteps()
{
  return(nLinearSteps);
}

const int& StrandBlockSolver::getNLevels()
{
  return(nLevels);
}

const int& StrandBlockSolver::getMgCycle()
{
  return(mgCycle);
}

const int& StrandBlockSolver::getNq()
{
  return(nq);
}

const int& StrandBlockSolver::getNFaces(){
  return(nFaces);
}

const int& StrandBlockSolver::getNPstr(){
  return(nPstr);
}

const int& StrandBlockSolver::getNBedges(){
  return(nBedges);
}

const int& StrandBlockSolver::getNDofs(){
  nDofs =(nFaces+nBedges)*(nPstr+2);
  return(nDofs);
}

const double& StrandBlockSolver::getConvLimit(){
  return(convLimit);
}

int* StrandBlockSolver::getFClip(){
  return(&fClip(0));
}

Array1D<int>* StrandBlockSolver::getFClipArray(){
  return(&fClip);
}

double* StrandBlockSolver::getQ(){
  return(&q(0,0,0));
}

Array3D<double>* StrandBlockSolver::getQArray(){
  return(&q);
}

double* StrandBlockSolver::getXc(){
  return(&xc(0,0,0));
}

double* StrandBlockSolver::getQx(){
  return(&qx(0,0,0,0));
}

double* StrandBlockSolver::getV(){
  return(&v(0,0));
}

Array2D<double>* StrandBlockSolver::getVArray(){
  return(&v);
}

double* StrandBlockSolver::getRms()
{
  return(&rms(0));
}

StrandSystem* StrandBlockSolver::getSystem()
{
  return(sys);
}

const int& StrandBlockSolver::getNGfaces(){
  return(nGfaces);
}

const int& StrandBlockSolver::getNFringe(){
  return(nFringe);
}

const int& StrandBlockSolver::getNPedges(){
  return(nPedges);
}

const int& StrandBlockSolver::getNEdges(){
  return(nEdges);
}

const int& StrandBlockSolver::getNqa()
{
  return(nqa);
}

const int& StrandBlockSolver::getPid()
{
  return(pid);
}

const int& StrandBlockSolver::getNdim()
{
  return(ndim);
}

const int& StrandBlockSolver::getNBpatches()
{
  return(nBpatches);
}

const int& StrandBlockSolver::getInviscid()
{
  return(inviscid);
}

const int& StrandBlockSolver::getViscous()
{
  return(viscous);
}

const int& StrandBlockSolver::getSource()
{
  return(source);
}

const int& StrandBlockSolver::getDissipation()
{
  return(dissipation);
}

int* StrandBlockSolver::getFTag(){
  return(&fTag(0));
}

Array1D<int>* StrandBlockSolver::getFTagArray(){
  return(&fTag);
}

int* StrandBlockSolver::getBTag(){
  return(&bTag(0));
}

Array1D<int>* StrandBlockSolver::getBTagArray(){
  return(&bTag);
}

int* StrandBlockSolver::getEdge(){
  return(&edge(0,0));
}

Array2D<int>* StrandBlockSolver::getEdgeArray(){
  return(&edge);
}

double* StrandBlockSolver::getFacu(){
  return(&facu(0,0,0));
}

Array3D<double>* StrandBlockSolver::getFacuArray(){
  return(&facu);
}

double* StrandBlockSolver::getFacs(){
  return(&facs(0,0,0));
}

Array3D<double>* StrandBlockSolver::getFacsArray(){
  return(&facs);
}

double* StrandBlockSolver::getXvu(){
  return(&xvu(0,0));
}

Array2D<double>* StrandBlockSolver::getXvuArray(){
  return(&xvu);
}

double* StrandBlockSolver::getXvs(){
  return(&xvs(0,0));
}

Array2D<double>* StrandBlockSolver::getXvsArray(){
  return(&xvs);
}

double* StrandBlockSolver::getR(){
  return(&r(0,0,0));
}

Array3D<double>* StrandBlockSolver::getRArray(){
  return(&r);
}

double* StrandBlockSolver::getQa(){
  return(&qa(0,0,0));
}

Array3D<double>* StrandBlockSolver::getQaArray(){
  return(&qa);
}

int* StrandBlockSolver::getLimFlag(){
  return(&limFlag);
}

int* StrandBlockSolver::getGradQFlag(){
  return(&gradQFlag);
}

int* StrandBlockSolver::getGradQaFlag(){
  return(&gradQaFlag);
}

int* StrandBlockSolver::getNodalQFlag(){
  return(&nodalQFlag);
}

int* StrandBlockSolver::getNodalQaFlag(){
  return(&nodalQaFlag);
}
