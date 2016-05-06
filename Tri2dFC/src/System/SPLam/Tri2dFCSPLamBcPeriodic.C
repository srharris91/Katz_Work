#include "Tri2dFCSPLamBcPeriodic.h"


// [Tri2dFCSPLamBcPeriodic]
Tri2dFCSPLamBcPeriodic::Tri2dFCSPLamBcPeriodic()
{
  cout << "\n***Periodic boundary conditions not yet implemented.***" << endl;
  exit(0);
}
// [Tri2dFCSPLamBcPeriodic]


// [~Tri2dFCSPLamBcPeriodic]
Tri2dFCSPLamBcPeriodic::~Tri2dFCSPLamBcPeriodic()
{
}
// [~Tri2dFCSPLamBcPeriodic]


// [BCVector]
void Tri2dFCSPLamBcPeriodic::BCVector(const double* nx,
				       const double* q,
				       const double* qa,
				       double* rb)
{
}
// [BCVector]


// [BCVectorSelfJacobian]
void Tri2dFCSPLamBcPeriodic::BCVectorSelfJacobian(const double* nx,
						   const double* q,
						   const double* qa,
						   double* M)
{
}
// [BCVectorSelfJacobian]


// [BCSelectionMatrix]
void Tri2dFCSPLamBcPeriodic::BCSelectionMatrix(const double* nx,
						const double* q,
						const double* qa,
						double* L)
{
}
// [BCSelectionMatrix]
