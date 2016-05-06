#include "Tri2dFCSPLamBcViscousWall.h"


// [Tri2dFCSPLamBcViscousWall]
Tri2dFCSPLamBcViscousWall::Tri2dFCSPLamBcViscousWall()
{
}
// [Tri2dFCSPLamBcViscousWall]


// [~Tri2dFCSPLamBcViscousWall]
Tri2dFCSPLamBcViscousWall::~Tri2dFCSPLamBcViscousWall()
{
}
// [~Tri2dFCSPLamBcViscousWall]


// [BCVector]
void Tri2dFCSPLamBcViscousWall::BCVector(const double* nx,
					 const double* q,
					 const double* qa,
					 double* rb)
{
  rb[0] = 0.             ;
  rb[1] = qa[1]          ; //x-velocity
  rb[2] = qa[2]          ; //y-velocity
  rb[3] = qa[3]-bValue[3]; //temperature
}
// [BCVector]


// [BCVectorSelfJacobian]
void Tri2dFCSPLamBcViscousWall::BCVectorSelfJacobian(const double* nx,
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

  M[0 ] = 0.;
  M[1 ] = 0.;
  M[2 ] = 0.;
  M[3 ] = 0.;

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
void Tri2dFCSPLamBcViscousWall::BCSelectionMatrix(const double* nx,
						const double* q,
						const double* qa,
						double* L)
{
  for (int n=0; n<nq*nq; n++) L[n] = 0.;
  L[0] = 1.; //mass
}
// [BCSelectionMatrix]


// [surfaceForces]
void Tri2dFCSPLamBcViscousWall::surfaceForces(const double* A,
					      const double* q,
					      const double* qx,
					      const double* qy,
					      const double* qa,
					      const double* qax,
					      const double* qay,
					      double* force)
{
  double mu,ux,vx,uy,vy,dd,sxx,syy,sxy,syx;
  mu       = qa [4];
  ux       = qax[1];
  vx       = qax[2];
  uy       = qay[1];
  vy       = qay[2];
  dd       =(ux+vy)/3.;
  sxx      = 2.*mu*(ux-dd);
  syy      = 2.*mu*(vy-dd);
  sxy      =    mu*(uy+vx);
  syx      = sxy;
  force[0] = qa[0]*A[0]-sxx*A[0]-sxy*A[1];
  force[1] = qa[0]*A[1]-syx*A[0]-syy*A[1];
}
// [surfaceForces]


// [surfaceSolution]
void Tri2dFCSPLamBcViscousWall::surfaceSolution(ofstream& ffile,
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
