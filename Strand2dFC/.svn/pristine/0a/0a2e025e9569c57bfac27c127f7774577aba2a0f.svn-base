#include "Strand2dFCSystem.h"


// [Strand2dFCSystem]
Strand2dFCSystem::Strand2dFCSystem()
{
  iPrint = 0;
  iTest = 0;
  iDebug = 0;
  iSolnFile = 0;
  iResdFile = 0;
  iErrFile = 0;
  nq = 0;
  nqa = 0;
  nqGradQ = 0;
  nqaGradQa = 0;
  ndim = 0;
  ncomp = 0;
  isolution = 0;
  nBpatches = 0;
  inviscid = 0;
  viscous = 0;
  source = 0;
  sourceMMS = 0;
  dissipation = 0;

}
// [Strand2dFCSystem]


// [~Strand2dFCSystem]
Strand2dFCSystem::~Strand2dFCSystem()
{
}
// [~Strand2dFCSystem]


// [rhsSource]
void Strand2dFCSystem::rhsSource(const int& npts,
				 const double* q,
				 const double* qa,
				 const double* qax,
				 const double* qay,
				 double* s)
{
  for (int n=0; n<npts*nq; n++) s[n] = 0.;
}
// [rhsSource]


// [lhsSourceJacobian]
void Strand2dFCSystem::lhsSourceJacobian(const int& npts,
					 const double* q,
					 const double* qa,
					 const double* qax,
					 const double* qay,
					 double* A)
{
  for (int n=0; n<npts*nq*nq; n++) A[n] = 0.;
}
// [lhsSourceJacobian]


// [initWallDist]
void Strand2dFCSystem::initWallDist(const int& npts,
                                    const double* dw,
                                    double* qa)
{  
}
// [initWallDist]
