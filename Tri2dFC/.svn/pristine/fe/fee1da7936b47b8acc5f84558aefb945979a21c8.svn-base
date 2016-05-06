#include "Tri2dFCSPLamBcInviscidWall.h"


// [Tri2dFCSPLamBcInviscidWall]
Tri2dFCSPLamBcInviscidWall::Tri2dFCSPLamBcInviscidWall()
{
}
// [Tri2dFCSPLamBcInviscidWall]


// [~Tri2dFCSPLamBcInviscidWall]
Tri2dFCSPLamBcInviscidWall::~Tri2dFCSPLamBcInviscidWall()
{
}
// [~Tri2dFCSPLamBcInviscidWall]


// [BCVector]
void Tri2dFCSPLamBcInviscidWall::BCVector(const double* nx,
					  const double* q,
					  const double* qa,
					  double* rb)
{
  double ru,rv,run;
  ru    = q[1];
  rv    = q[2];
  run   = nx[0]*ru+nx[1]*rv;

  rb[0] = 0.;
  rb[1] = 0.;
  rb[2] = 0.;
  rb[3] = run-0.; //normal velocity
  //rb[3] = 0.;
}
// [BCVector]


// [BCVectorSelfJacobian]
void Tri2dFCSPLamBcInviscidWall::BCVectorSelfJacobian(const double* nx,
						      const double* q,
						      const double* qa,
						      double* M)
{
  double rr,u,v,un;
  rr    = 1./q[0];
  u     = qa[1];
  v     = qa[2];
  un    = nx[0]*u+nx[1]*v;

  M[0 ] = 0.;
  M[1 ] = 0.;
  M[2 ] = 0.;
  M[3 ] = 0.;

  M[4 ] = 0.;
  M[5 ] = 0.;
  M[6 ] = 0.;
  M[7 ] = 0.;

  M[8 ] = 0.;
  M[9 ] = 0.;
  M[10] = 0.;
  M[11] = 0.;

  M[12] = 0.;
  M[13] = nx[0];
  M[14] = nx[1];
  M[15] = 0.;

  /*
  M[12] = 0.;
  M[13] = 0.;
  M[14] = 0.;
  M[15] = 0.;
  */
}
// [BCVectorSelfJacobian]


// [BCSelectionMatrix]
void Tri2dFCSPLamBcInviscidWall::BCSelectionMatrix(const double* nx,
						   const double* q,
						   const double* qa,
						   double* L)
{
  double ccr,qq,u,v,ut,un,c;
  u     = qa[1];
  v     = qa[2];
  un    = nx[0]*u+nx[1]*v;
  ut    =-nx[1]*u+nx[0]*v;
  qq    = u*u+v*v;
  c     = sqrt(gamma*qa[0]/q[0]);
  ccr   = 1./(c*c);

  L[0 ] =-ccr*gm1*qq*.5+1.;
  L[1 ] = ccr*gm1*u;
  L[2 ] = ccr*gm1*v;
  L[3 ] =-ccr*gm1;

  L[4 ] = ut;
  L[5 ] = nx[1];
  L[6 ] =-nx[0];
  L[7 ] = 0.;

  L[8 ] = ccr*( gm1*qq*.5-c*un   )*.5;
  L[9 ] = ccr*(-gm1*u    +c*nx[0])*.5;
  L[10] = ccr*(-gm1*v    +c*nx[1])*.5;
  L[11] = ccr*( gm1              )*.5;

  L[12] = 0.;
  L[13] = 0.;
  L[14] = 0.;
  L[15] = 0.;

  /*//solves mass, tangential momentum, and energy
  L[0 ] = 1.;
  L[1 ] = 0.;
  L[2 ] = 0.;
  L[3 ] = 0.;

  L[4 ] = 0.;
  L[5 ] =-nx[1];
  L[6 ] = nx[0];
  L[7 ] = 0.;

  L[8 ] = 0.;
  L[9 ] = 0.;
  L[10] = 0.;
  L[11] = 1.;

  L[12] = 0.;
  L[13] = 0.;
  L[14] = 0.;
  L[15] = 0.;
  */

  /*//solves all interior equations
  L[0 ] = 1.;
  L[1 ] = 0.;
  L[2 ] = 0.;
  L[3 ] = 0.;

  L[4 ] = 0.;
  L[5 ] = 1.;
  L[6 ] = 0.;
  L[7 ] = 0.;

  L[8 ] = 0.;
  L[9 ] = 0.;
  L[10] = 1.;
  L[11] = 0.;

  L[12] = 0.;
  L[13] = 0.;
  L[14] = 0.;
  L[15] = 1.;
  */
}
// [BCSelectionMatrix]


// [surfaceForces]
void Tri2dFCSPLamBcInviscidWall::surfaceForces(const double* A,
					       const double* q,
					       const double* qx,
					       const double* qy,
					       const double* qa,
					       const double* qax,
					       const double* qay,
					       double* force)
{
  force[0] = qa[0]*A[0];
  force[1] = qa[0]*A[1];
}
// [surfaceForces]


// [surfaceSolution]
void Tri2dFCSPLamBcInviscidWall::surfaceSolution(ofstream& ffile,
						 const double* x,
						 const double* y,
						 const double* q,
						 const double* qa,
						 const double* qx,
						 const double* qy,
						 const double* qax,
						 const double* qay)
{
  double p0,u0,v0,t0,r0,q0,c0,Cp;
  p0 = bValue[0];
  u0 = bValue[1];
  v0 = bValue[2];
  t0 = bValue[3];
  r0 = p0/(rGas*t0);
  q0 = u0*u0+v0*v0;
  c0 = 1.; //assume chord length of 1 for now.
  Cp =(qa[0]-p0)/(.5*r0*q0*c0);
  ffile << x[0] << " " << -Cp << endl;
}
// [surfaceSolution]
