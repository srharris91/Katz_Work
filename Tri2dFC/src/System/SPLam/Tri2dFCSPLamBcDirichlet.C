#include "Tri2dFCSPLamBcDirichlet.h"


// [Tri2dFCSPLamBcDirichlet]
Tri2dFCSPLamBcDirichlet::Tri2dFCSPLamBcDirichlet()
{
}
// [Tri2dFCSPLamBcDirichlet]


// [~Tri2dFCSPLamBcDirichlet]
Tri2dFCSPLamBcDirichlet::~Tri2dFCSPLamBcDirichlet()
{
}
// [~Tri2dFCSPLamBcDirichlet]


// [BCVector]
void Tri2dFCSPLamBcDirichlet::BCVector(const double* nx,
				       const double* q,
				       const double* qa,
				       double* rb)
{
  rb[0] = qa[0]-bValue[0]; //pressure
  rb[1] = qa[1]-bValue[1]; //x-velocity
  rb[2] = qa[2]-bValue[2]; //y-velocity
  rb[3] = qa[3]-bValue[3]; //temperature
}
// [BCVector]


// [BCVectorSelfJacobian]
void Tri2dFCSPLamBcDirichlet::BCVectorSelfJacobian(const double* nx,
						   const double* q,
						   const double* qa,
						   double* M)
{
  double rr,e,u,v,qq,grr;
  rr    = 1./q[0];
  e     = q[3]*rr;
  u     = qa[1];
  v     = qa[2];
  qq    = u*u+v*v;
  grr   = gm1*rr/rGas;

  M[0 ] = gm1*qq*.5;
  M[1 ] =-gm1*u;
  M[2 ] =-gm1*v;
  M[3 ] = gm1;

  M[4 ] =-rr*u;
  M[5 ] = rr;
  M[6 ] = 0.;
  M[7 ] = 0.;

  M[8 ] =-rr*v;
  M[9 ] = 0.;
  M[10] = rr;
  M[11] = 0.;

  M[12] = grr*(qq-e);
  M[13] =-grr*u;
  M[14] =-grr*v;
  M[15] = grr;
}
// [BCVectorSelfJacobian]


// [BCSelectionMatrix]
void Tri2dFCSPLamBcDirichlet::BCSelectionMatrix(const double* nx,
						const double* q,
						const double* qa,
						double* L)
{
  for (int n=0; n<nq*nq; n++) L[n] = 0.; //none
}
// [BCSelectionMatrix]
