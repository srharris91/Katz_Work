#include "StrandSPTurbSABcDirichlet.h"


// [StrandSPTurbSABcDirichlet]
StrandSPTurbSABcDirichlet::StrandSPTurbSABcDirichlet()
{
}
// [StrandSPTurbSABcDirichlet]


// [~StrandSPTurbSABcDirichlet]
StrandSPTurbSABcDirichlet::~StrandSPTurbSABcDirichlet()
{
}
// [~StrandSPTurbSABcDirichlet]


// [BCVector]
void StrandSPTurbSABcDirichlet::BCVector(const double* nx,
					 const double* wx,
					 const double* qe,
					 const double* qae,
					 const double* q,
					 const double* qa,
					 double* r)
{
  r[0] = qa[0]-bValue[0];
  r[1] = qa[1]-bValue[1];
  r[2] = qa[2]-bValue[2];
  r[3] = qa[3]-bValue[3];
  r[4] = qa[4]-bValue[4];
}
// [BCVector]


// [BCVectorSelfJacobian]
void StrandSPTurbSABcDirichlet::BCVectorSelfJacobian(const double* nx,
						     const double* qe,
						     const double* qae,
						     const double* q,
						     const double* qa,
						     double* M)
{
  double rr,e,u,v,qq,grr,nu;
  rr    = 1./q[0];
  e     = q[3]*rr;
  u     = qa[1];
  v     = qa[2];
  nu    = qa[4];
  qq    = u*u+v*v;
  grr   = gm1*rr/rGas;

  M[0 ] = gm1*qq*.5;
  M[1 ] =-gm1*u;
  M[2 ] =-gm1*v;
  M[3 ] = gm1;
  M[4 ] = 0.;

  M[5 ] =-rr*u;
  M[6 ] = rr;
  M[7 ] = 0.;
  M[8 ] = 0.;
  M[9 ] = 0.;

  M[10] =-rr*v;
  M[11] = 0.;
  M[12] = rr;
  M[13] = 0.;
  M[14] = 0.;

  M[15] = grr*(qq-e);
  M[16] =-grr*u;
  M[17] =-grr*v;
  M[18] = grr;
  M[19] = 0.;
  
  M[20] =-rr*nu;
  M[21] = 0.;
  M[22] = 0.;
  M[23] = 0.;
  M[24] = rr;
  
}
// [BCVectorSelfJacobian]


// [BCVectorInteriorJacobian]
void StrandSPTurbSABcDirichlet::BCVectorInteriorJacobian(const double* nx,
						      const double* qe,
						      const double* qae,
						      const double* q,
						      const double* qa,
						      double* M)
{
  for (int n=0; n<nq*nq; n++) M[n] = 0.;
}
// [BCVectorInteriorJacobian]
