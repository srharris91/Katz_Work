#include "StrandSPLamBcViscousWall.h"


// [StrandSPLamBcViscousWall]
StrandSPLamBcViscousWall::StrandSPLamBcViscousWall()
{
}
// [StrandSPLamBcViscousWall]


// [~StrandSPLamBcViscousWall]
StrandSPLamBcViscousWall::~StrandSPLamBcViscousWall()
{
}
// [~StrandSPLamBcViscousWall]


// [surfaceForces]
void StrandSPLamBcViscousWall::surfaceForces(const double* A,
					     const double* q,
					     const double* qa,
					     const double* qx,
					     const double* qax,
					     double* force)
{
  double mu,ux,vx,uy,vy,dd,sxx,syy,sxy,syx;

  mu       = qa [4];
  ux       = qax[1];
  vx       = qax[2];
  uy       = qax[nqa+1];
  vy       = qax[nqa+2];
  dd       =(ux+vy)/3.;
  sxx      = 2.*mu*(ux-dd);
  syy      = 2.*mu*(vy-dd);
  sxy      =    mu*(uy+vx);
  syx      = sxy;
  force[0] =-qa[0]*A[0]+sxx*A[0]+sxy*A[1];
  force[1] =-qa[0]*A[1]+syx*A[0]+syy*A[1];
}
// [surfaceForces]


// [BCVector]
void StrandSPLamBcViscousWall::BCVector(const double* nx,
					const double* wx,
					const double* qe,
					const double* qae,
					const double* q,
					const double* qa,
					double* r)
{
  double pb,ub,vb,tb,pe,ue,ve,te;

  pb   = qa[0];
  ub   = qa[1];
  vb   = qa[2];
  tb   = qa[3];

  pe   = qae[0];
  ue   = qae[1];
  ve   = qae[2];
  te   = qae[3];

  r[0] = pb-pe;
  r[1] = ub+ue-2.*wx[0];
  r[2] = vb+ve-2.*wx[1];
  r[3] = tb-te;
}
// [BCVector]


// [BCVectorSelfJacobian]
void StrandSPLamBcViscousWall::BCVectorSelfJacobian(const double* nx,
						    const double* qe,
						    const double* qae,
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


// [BCVectorInteriorJacobian]
void StrandSPLamBcViscousWall::BCVectorInteriorJacobian(const double* nx,
							const double* qe,
							const double* qae,
							const double* q,
							const double* qa,
							double* M)
{
  double rr,e,u,v,qq,grr;
  rr    = 1./qe[0];
  e     = qe[3]*rr;
  u     = qae[1];
  v     = qae[2];
  qq    = u*u+v*v;
  grr   = gm1*rr/rGas;

  M[0 ] =-gm1*qq*.5;
  M[1 ] = gm1*u;
  M[2 ] = gm1*v;
  M[3 ] =-gm1;

  M[4 ] =-rr*u;
  M[5 ] = rr;
  M[6 ] = 0.;
  M[7 ] = 0.;

  M[8 ] =-rr*v;
  M[9 ] = 0.;
  M[10] = rr;
  M[11] = 0.;

  M[12] =-grr*(qq-e);
  M[13] = grr*u;
  M[14] = grr*v;
  M[15] =-grr;
}
// [BCVectorInteriorJacobian]
