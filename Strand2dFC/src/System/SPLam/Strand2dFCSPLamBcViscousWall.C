#include "Strand2dFCSPLamBcViscousWall.h"


// [Strand2dFCSPLamBcViscousWall]
Strand2dFCSPLamBcViscousWall::Strand2dFCSPLamBcViscousWall()
{
}
// [Strand2dFCSPLamBcViscousWall]


// [~Strand2dFCSPLamBcViscousWall]
Strand2dFCSPLamBcViscousWall::~Strand2dFCSPLamBcViscousWall()
{
}
// [~Strand2dFCSPLamBcViscousWall]


// [BCVector]
void Strand2dFCSPLamBcViscousWall::BCVector(const double* nx,
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
void Strand2dFCSPLamBcViscousWall::BCVectorSelfJacobian(const double* nx,
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
void Strand2dFCSPLamBcViscousWall::BCSelectionMatrix(const double* nx,
						     const double* q,
						     const double* qa,
						     double* L)
{
  for (int n=0; n<nq*nq; n++) L[n] = 0.;
  L[0] = 1.; //mass
}
// [BCSelectionMatrix]


// [BCPenalty]
void Strand2dFCSPLamBcViscousWall::BCPenalty(const int& inout,
					     const double* A,
					     const double& Pinv0,
					     const double* q,
					     const double* qa,
					     const double* g,
					     const double* gv,
					     const double* uw,
					     double* rb)
{
  double Ax,Ay,ds,Nx,Ny,Tx,Ty,rl,rul,rvl,re,rqq,p,rhl,rr,rur,rvr,rhr,dd,mu,
    u,v,h,qq,cc,ccr,ccr2,c,ut,un,m,l[nq],R[nq*nq],S[nq*nq],M[nq*nq],qw[nq];

  Ax      = A[0];
  Ay      = A[1];
  ds      = sqrt(Ax*Ax+Ay*Ay);
  Nx      = Ax/ds;
  Ny      = Ay/ds;
  Tx      = Ny;
  Ty      =-Nx;

  rl      = q[0];
  rul     = q[1];
  rvl     = q[2];
  re      = q[3];
  rqq     =(rul*rul+rvl*rvl)/rl;
  p       = gm1*(re-.5*rqq);
  rhl     = re+p;

  // form the inviscid wall state by subtracting off normal velocity
  ut      = rul*Tx+rvl*Ty;
  qw[0]   = q[0];
  qw[1]   = Ny*ut;
  qw[2]   =-Nx*ut;
  qw[3]   = q[3];
  //qw[3]   = rGas/gm1*q[0]*bValue[3]+.5*q[0]*(ut*ut); //temperature condition

  rr      = qw[0];
  rur     = qw[1];
  rvr     = qw[2];
  re      = qw[3];
  rqq     =(rur*rur+rvr*rvr)/rr;
  p       = gm1*(re-.5*rqq);
  rhr     = re+p;

  rl      = sqrt(rl);
  rr      = sqrt(rr);
  dd      = 1./(rl+rr);
  rl      = 1./rl;
  rr      = 1./rr;
  u       =(rul*rl+rur*rr)*dd;
  v       =(rvl*rl+rvr*rr)*dd;
  h       =(rhl*rl+rhr*rr)*dd;
  qq      = .5*(u*u+v*v);
  cc      = gm1*(h-qq);
  ccr     = 1./cc;
  ccr2    = .5*ccr;
  c       = sqrt(cc);
  ut      = u*Tx+v*Ty;
  un      = u*Nx+v*Ny;

  if (inout > 0){
    l[0]    = 0.;
    l[1]    = 0.;
    l[2]    = ds*c;
    l[3]    = 0.;
  }
  else{
    l[0]    = 0.;
    l[1]    = 0.;
    l[2]    = 0.;
    l[3]    =-ds*c;
  }

  R[0 ] = l[0];
  R[1 ] = 0.;
  R[2 ] = l[2];
  R[3 ] = l[3];
  
  R[4 ] = l[0]*u;
  R[5 ] = l[1]*Tx;
  R[6 ] = l[2]*(u+Nx*c);
  R[7 ] = l[3]*(u-Nx*c);
  
  R[8 ] = l[0]*v;
  R[9 ] = l[1]*Ty;
  R[10] = l[2]*(v+Ny*c);
  R[11] = l[3]*(v-Ny*c);
  
  R[12] = l[0]*qq;
  R[13] = l[1]*ut;
  R[14] = l[2]*(h+un*c);
  R[15] = l[3]*(h-un*c);
  
  S[0 ] =-gm1*ccr*qq+1.;
  S[1 ] = gm1*ccr*u;
  S[2 ] = gm1*ccr*v;
  S[3 ] =-gm1*ccr;
  
  S[4 ] =-ut;
  S[5 ] = Tx;
  S[6 ] = Ty;
  S[7 ] = 0.;
  
  S[8 ] = ccr2*(gm1*qq-c*un);
  S[9 ] =-ccr2*(gm1*u -c*Nx);
  S[10] =-ccr2*(gm1*v -c*Ny);
  S[11] = ccr2* gm1;
  
  S[12] = ccr2*(gm1*qq+c*un);
  S[13] =-ccr2*(gm1*u +c*Nx);
  S[14] =-ccr2*(gm1*v +c*Ny);
  S[15] = ccr2* gm1;
  
  for (int k=0; k<nq; k++) l[k] = qw[k]-q[k];
  
  matmul(nq,nq,nq,&R[0],&S[0],&M[0] );
  matmul(nq,nq,1 ,&M[0],&l[0],&rb[0]);

  for (int k=0; k<nq; k++) rb[k] *= 2.; //factor of two for well-posedness

  // add in viscous penalty
  qw[0] = q[0];
  qw[1] = 0.;
  qw[2] = 0.;
  qw[3] = rGas/gm1*q[0]*bValue[3]; //temperature condition
  rl    = q[0];
  mu    = qa[4];
  p     = .25*Pinv0*max(gamma*mu/(rl*Prn),5.*mu/(3.*rl));
  for (int k=0; k<nq; k++) rb[k] += p*(qw[k]-q[k]);

  //if (viscous) for (int k=0; k<nq; k++) rb[k] -= gv[k];
}
// [BCPenalty]


// [surfaceForces]
void Strand2dFCSPLamBcViscousWall::surfaceForces(const double* xs,
						 const double* ys,
						 const double* q,
						 const double* qa,
						 const double* qx,
						 const double* qy,
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
  force[0] = qa[0]*ys[0]-sxx*ys[0]+sxy*xs[0];
  force[1] =-qa[0]*xs[0]-syx*ys[0]+syy*xs[0];
}
// [surfaceForces]
