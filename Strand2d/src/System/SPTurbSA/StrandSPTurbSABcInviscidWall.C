#include "StrandSPTurbSABcInviscidWall.h"


// [StrandSPTurbSABcInviscidWall]
StrandSPTurbSABcInviscidWall::StrandSPTurbSABcInviscidWall()
{
}
// [StrandSPTurbSABcInviscidWall]


// [~StrandSPTurbSABcInviscidWall]
StrandSPTurbSABcInviscidWall::~StrandSPTurbSABcInviscidWall()
{
}
// [~StrandSPTurbSABcInviscidWall]


// [surfaceForces]
void StrandSPTurbSABcInviscidWall::surfaceForces(const double* A,
						 const double* q,
						 const double* qa,
						 const double* qx,
						 const double* qax,
						 double* force)
{
  force[0] =-qa[0]*A[0];
  force[1] =-qa[0]*A[1];
}
// [surfaceForces]


// [BCVector]
void StrandSPTurbSABcInviscidWall::BCVector(const double* nx,
					    const double* wx,
					    const double* qe,
					    const double* qae,
					    const double* q,
					    const double* qa,
					    double* r)
{
  double Nx,Ny,Tx,Ty,ute,une,utb,unb,unw,du,dv;

  Nx   = nx[0];
  Ny   = nx[1];
  Tx   = nx[1];
  Ty   =-nx[0];

  ute  = Tx*qae[1]+Ty*qae[2];
  une  = Nx*qae[1]+Ny*qae[2];
  utb  = Tx*qa [1]+Ty*qa [2];
  unb  = Nx*qa [1]+Ny*qa [2];
  unw  = Nx*wx [0]+Ny*wx [1];

  //du   = Tx*(utb-ute)+Nx*(unb+une-2.*unw);
  //dv   = Ty*(utb-ute)+Ny*(unb+une-2.*unw);

  du   = Ny*(utb-ute)-Ty*(unb+une-2.*unw);
  dv   =-Nx*(utb-ute)+Tx*(unb+une-2.*unw);

  r[0] = q[0]-qe[0]; //density
  r[1] = du;         //tangential velocity
  r[2] = dv;         //normal velocity
  r[3] = q[3]-qe[3]; //total energy per unit volume
  r[4] = q[4]-qe[4]; //turbulent eddy viscosity
}
// [BCVector]


// [BCVectorSelfJacobian]
void StrandSPTurbSABcInviscidWall::BCVectorSelfJacobian(const double* nx,
						     const double* qe,
						     const double* qae,
						     const double* q,
						     const double* qa,
						     double* M)
{
  double rr,u,v;

  rr    = 1./q[0];
  u     = qa[1];
  v     = qa[2];

  M[0 ] = 1.;
  M[1 ] = 0.;
  M[2 ] = 0.;
  M[3 ] = 0.;
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

  M[15] = 0.;
  M[16] = 0.;
  M[17] = 0.;
  M[18] = 1.;
  M[19] = 0.;
  
  M[20] = 0.;
  M[21] = 0.;
  M[22] = 0.;
  M[23] = 0.;
  M[24] = 1.;
}
// [BCVectorSelfJacobian]


// [BCVectorInteriorJacobian]
void StrandSPTurbSABcInviscidWall::BCVectorInteriorJacobian(const double* nx,
							 const double* qe,
							 const double* qae,
							 const double* q,
							 const double* qa,
							 double* M)
{
  double Nx,Ny,Tx,Ty,rr,ut,un,u,v,a,b,c;

  Nx   = nx[0];
  Ny   = nx[1];
  Tx   = nx[1];
  Ty   =-nx[0];
  rr   = 1./qe[0];
  u    = qae[1];
  v    = qae[2];
  ut   = Tx*u+Ty*v;
  un   = Nx*u+Ny*v;
  a    = Nx*Ty+Ny*Tx;
  b    = 2.*Nx*Tx;
  c    = 2.*Ny*Ty;

  M[0 ] =-1.;
  M[1 ] = 0.;
  M[2 ] = 0.;
  M[3 ] = 0.;
  M[4 ] = 0.;

  M[5 ] = rr*(a*u+c*v);
  M[6 ] =-rr*a;
  M[7 ] =-rr*c;
  M[8 ] = 0.;
  M[9 ] = 0.;


  M[10] =-rr*(b*u+a*v);
  M[11] = rr*b;
  M[12] = rr*a;
  M[13] = 0.;
  M[14] = 0.;

  M[15] = 0.;
  M[16] = 0.;
  M[17] = 0.;
  M[18] =-1.;
  M[19] = 0.;
  
  M[20] = 0.;
  M[21] = 0.;
  M[22] = 0.;
  M[23] = 0.;
  M[24] =-1.;
}
// [BCVectorInteriorJacobian]
