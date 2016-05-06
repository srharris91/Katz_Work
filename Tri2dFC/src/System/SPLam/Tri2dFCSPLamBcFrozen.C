#include "Tri2dFCSPLamBcFrozen.h"


// [Tri2dFCSPLamBcFrozen]
Tri2dFCSPLamBcFrozen::Tri2dFCSPLamBcFrozen()
{
}
// [Tri2dFCSPLamBcFrozen]


// [~Tri2dFCSPLamBcFrozen]
Tri2dFCSPLamBcFrozen::~Tri2dFCSPLamBcFrozen()
{
}
// [~Tri2dFCSPLamBcFrozen]


// [BCVector]
void Tri2dFCSPLamBcFrozen::BCVector(const double* nx,
				    const double* q,
				    const double* qa,
				    double* rb)
{
  for (int n=0; n<nq; n++) rb[n] = 0.; //none
}
// [BCVector]


// [BCVectorSelfJacobian]
void Tri2dFCSPLamBcFrozen::BCVectorSelfJacobian(const double* nx,
						const double* q,
						const double* qa,
						double* M)
{
  for (int n=0; n<nq*nq; n++) M[n     ] = 0.;
  for (int n=0; n<nq   ; n++) M[n*nq+n] = 1.; //identity
}
// [BCVectorSelfJacobian]


// [BCSelectionMatrix]
void Tri2dFCSPLamBcFrozen::BCSelectionMatrix(const double* nx,
					     const double* q,
					     const double* qa,
					     double* L)
{
  for (int n=0; n<nq*nq; n++) L[n] = 0.; //none
}
// [BCSelectionMatrix]
