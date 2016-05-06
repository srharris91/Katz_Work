#include "StrandSystem.h"


// [StrandSystem]
StrandSystem::StrandSystem()
{
  iPrint = 0;
  iTest = 0;
  iDebug = 0;
  nq = 0;
  nqa = 0;
  ndim = 0;
  ncomp = 0;
  nBpatches = 0;
  isolution = 0;
  inviscid = 0;
  viscous = 0;
  source = 0;
  dissipation = 0;
}
// [StrandSystem]


// [~StrandSystem]
StrandSystem::~StrandSystem()
{
}
// [~StrandSystem]


// [inputRead]
void StrandSystem::inputRead(const string& inputFile)
{
  cout << "\ninputRead not redefined for StrandSystem sub-class. Terminating."
       << endl;
  exit(0);
}
// [inputRead]


// [prepSetup]
void StrandSystem::prepSetup(const int& iPrint,
			     const int& iTest,
			     const int& iDebug,
			     const int& tmp,
			     int& nq,
			     int& nqa,
			     int& ndim,
			     int& inviscid,
			     int& viscous,
			     int& source,
			     int& sourceMMS,
			     int& dissipation,
			     int& nBpatches,
			     int* iqgradT,
			     int* iqagradT,
			     double* dlim,
			     double* rmsNorm)
{
  cout << "\nprepSetup not redefined for StrandSystem sub-class. Terminating."
       << endl;
  exit(0);
}
// [prepSetup]


// [initWallDist]
void StrandSystem::initWallDist(const int& npts,
				const double* dw,
				double* qa)
{
}
// [initWallDist]


// [rhsSource]
void StrandSystem::rhsSource(const int& npts,
			     const double* q,
			     const double* qa,
			     const double* qx,
			     const double* qax,
			     double* f)
{
  for (int n=0; n<nq*npts; n++) f[n] = 0.;
}
// [rhsSource]


// [rhsSourceCoarse]
void StrandSystem::rhsSourceCoarse(const int& npts,
				   const double* v,
				   const double* q,
				   const double* qa,
				   double* f)
{
  for (int n=0; n<nq*npts; n++) f[n] = 0.;
}
// [rhsSourceCoarse]


// [lhsVisFluxJacobian]
void StrandSystem::lhsVisFluxJacobian(const int& npts,
				      const double* A,
				      const double* q,
				      const double* qa,
				      const double* qx,
				      const double* qy,
				      const double* qax,
				      const double* qay,
				      double* M)
{
  for (int n=0; n<npts*nq*nq; n++) M[n] = 0.;
}
// [lhsVisFluxJacobian]


// [lhsSourceJacobian]
void StrandSystem::lhsSourceJacobian(const int& npts,
				     const double* v,
				     const double* q,
				     const double* qa,
				     const double* qx,
				     const double* qax,
				     double* M)
{
  for (int n=0; n<npts*nq*nq; n++) M[n] = 0.;
}
// [lhsSourceJacobian]


// [lhsSourceJacobianCoarse]
void StrandSystem::lhsSourceJacobianCoarse(const int& npts,
					   const double* v,
					   const double* q,
					   const double* qa,
					   double* M)
{
  for (int n=0; n<nq*nq*npts; n++) M[n] = 0.;
}
// [rhsSourceJacobianCoarse]


// [lhsPreconJacobian]
void StrandSystem::lhsPreconJacobian(const int& npts,
				     const double* q,
				     const double* qa,
				     double* M)
{
  int iM;
  for (int n=0; n<npts*nq*nq; n++) M[n] = 0.;
  for (int n=0; n<npts; n++){
    iM = nq*nq*n;
    for (int i=0; i<nq; i++) M[iM+nq*i+i] = 1.;
  }
}
// [lhsPreconJacobian]


// [lhsConsVarJacobian]
void StrandSystem::lhsConsVarJacobian(const int& npts,
				      const double* q,
				      const double* qa,
				      double* M)
{
  int iM;
  for (int n=0; n<npts*nq*nq; n++) M[n] = 0.;
  for (int n=0; n<npts; n++){
    iM = nq*nq*n;
    for (int i=0; i<nq; i++) M[iM+nq*i+i] = 1.;
  }
}
// [lhsConsVarJacobian]


// [lhsPrimVarJacobian]
void StrandSystem::lhsPrimVarJacobian(const int& npts,
				      const double* q,
				      const double* qa,
				      double* M)
{
  int iM;
  for (int n=0; n<npts*nq*nq; n++) M[n] = 0.;
  for (int n=0; n<npts; n++){
    iM = nq*nq*n;
    for (int i=0; i<nq; i++) M[iM+nq*i+i] = 1.;
  }
}
// [lhsPrimVarJacobian]
