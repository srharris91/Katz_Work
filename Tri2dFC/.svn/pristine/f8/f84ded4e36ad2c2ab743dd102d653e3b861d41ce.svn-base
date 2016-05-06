#include "Tri2dFCSPLamBcNothing.h"


// [Tri2dFCSPLamBcNothing]
Tri2dFCSPLamBcNothing::Tri2dFCSPLamBcNothing()
{
}
// [Tri2dFCSPLamBcNothing]


// [~Tri2dFCSPLamBcNothing]
Tri2dFCSPLamBcNothing::~Tri2dFCSPLamBcNothing()
{
}
// [~Tri2dFCSPLamBcNothing]


// [BCVector]
void Tri2dFCSPLamBcNothing::BCVector(const double* nx,
				       const double* q,
				       const double* qa,
				       double* rb)
{
  for (int n=0; n<nq; n++) rb[n] = 0.; //none
}
// [BCVector]


// [BCVectorSelfJacobian]
void Tri2dFCSPLamBcNothing::BCVectorSelfJacobian(const double* nx,
						   const double* q,
						   const double* qa,
						   double* M)
{
  for (int n=0; n<nq*nq; n++) M[n] = 0.;
}
// [BCVectorSelfJacobian]


// [BCSelectionMatrix]
void Tri2dFCSPLamBcNothing::BCSelectionMatrix(const double* nx,
						const double* q,
						const double* qa,
						double* L)
{
  for (int n=0; n<nq*nq; n++) L[n     ] = 0.;
  for (int n=0; n<nq   ; n++) L[n*nq+n] = 1.; //all EOMs
}
// [BCSelectionMatrix]
