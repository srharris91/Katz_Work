#include "Strand2dFCSPTurbSABcViscousWall.h"


// [Strand2dFCSPTurbSABcViscousWall]
Strand2dFCSPTurbSABcViscousWall::Strand2dFCSPTurbSABcViscousWall()
{
}
// [Strand2dFCSPTurbSABcViscousWall]


// [~Strand2dFCSPTurbSABcViscousWall]
Strand2dFCSPTurbSABcViscousWall::~Strand2dFCSPTurbSABcViscousWall()
{
}
// [~Strand2dFCSPTurbSABcViscousWall]


// [BCVector]
void Strand2dFCSPTurbSABcViscousWall::BCVector(const double* nx,
					       const double* q,
					       const double* qa,
					       double* rb)
{
  rb[0] = 0.             ;
  rb[1] = qa[1]          ; //x-velocity
  rb[2] = qa[2]          ; //y-velocity
  rb[3] = qa[3]-bValue[3]; //temperature
  rb[4] = qa[4]          ; //nu_tilde


}
// [BCVector]




// [BCVectorSelfJacobian]
void Strand2dFCSPTurbSABcViscousWall::BCVectorSelfJacobian(const double* nx,
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

  M[0 ] = 0.;
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


// [BCSelectionMatrix]
void Strand2dFCSPTurbSABcViscousWall::BCSelectionMatrix(const double* nx,
						        const double* q,
						        const double* qa,
						        double* L)
{
  for (int n=0; n<nq*nq; n++) L[n] = 0.;
  L[0] = 1.; //mass
}
// [BCSelectionMatrix]


// [BCPenalty]
void Strand2dFCSPTurbSABcViscousWall::BCPenalty(const int& inout,
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
    u,v,h,qq,cc,ccr,ccr2,c,ut,un,m,l[nq],R[nq*nq],S[nq*nq],M[nq*nq],qw[nq],
    rnur,rnul,nu;

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
  rnul    = q[4];
  rqq     =(rul*rul+rvl*rvl)/rl;
  p       = gm1*(re-.5*rqq);
  rhl     = re+p;

  // form the inviscid wall state by subtracting off normal velocity
  ut      = rul*Tx+rvl*Ty;
  qw[0]   = q[0];
  qw[1]   = Ny*ut;
  qw[2]   =-Nx*ut;
  qw[3]   = q[3];
  qw[4]   = q[4];
  //qw[3]   = rGas/gm1*q[0]*bValue[3]+.5*q[0]*(ut*ut); //temperature condition

  rr      = qw[0];
  rur     = qw[1];
  rvr     = qw[2];
  re      = qw[3];
  rnur    = qw[4];
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
  nu      =(rnul*rl+rnur*rr)*dd;
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
    l[4]    = 0.;
  }
  else{
    l[0]    = 0.;
    l[1]    = 0.;
    l[2]    = 0.;
    l[3]    =-ds*c;
    l[4]    = 0.;
  }

  R[0 ] = l[0];
  R[1 ] = 0.;
  R[2 ] = l[2];
  R[3 ] = l[3];
  R[4 ] = 0.;

  R[5 ] = l[0]*u;
  R[6 ] = l[1]*Tx;
  R[7 ] = l[2]*(u+Nx*c);
  R[8 ] = l[3]*(u-Nx*c);
  R[9 ] = 0.;

  R[10] = l[0]*v;
  R[11] = l[1]*Ty;
  R[12] = l[2]*(v+Ny*c);
  R[13] = l[3]*(v-Ny*c);
  R[14] = 0.;
   
  R[15] = l[0]*qq;
  R[16] = l[1]*ut;
  R[17] = l[2]*(h+un*c);
  R[18] = l[3]*(h-un*c);
  R[19] = 0.;

  R[20] = 0.;
  R[21] = 0.;
  R[22] = l[2]*nu;
  R[23] = l[3]*nu;
  R[24] = l[4];

  S[0 ] =-gm1*ccr*qq+1.;
  S[1 ] = gm1*ccr*u;
  S[2 ] = gm1*ccr*v;
  S[3 ] =-gm1*ccr;
  S[4 ] = 0.;
  
  S[5 ] =-ut;
  S[6 ] = Tx;
  S[7 ] = Ty;
  S[8 ] = 0.;
  S[9 ] = 0.;

  S[10] = ccr2*(gm1*qq-c*un);
  S[11] =-ccr2*(gm1*u -c*Nx);
  S[12] =-ccr2*(gm1*v -c*Ny);
  S[13] = ccr2* gm1;
  S[14] = 0.;

  S[15] = ccr2*(gm1*qq+c*un);
  S[16] =-ccr2*(gm1*u +c*Nx);
  S[17] =-ccr2*(gm1*v +c*Ny);
  S[18] = ccr2* gm1;
  S[19] = 0.;
  
  S[20] =-nu*gm1*ccr*qq;
  S[21] = nu*gm1*ccr*u;
  S[22] = nu*gm1*ccr*v;
  S[23] =-nu*gm1*ccr;
  S[24] = 1.;

  for (int k=0; k<nq; k++) l[k] = qw[k]-q[k];
  
  matmul(nq,nq,nq,&R[0],&S[0],&M[0] );
  matmul(nq,nq,1 ,&M[0],&l[0],&rb[0]);

  for (int k=0; k<nq; k++) rb[k] *= 2.; //factor of two for well-posedness

  // add in viscous penalty
  qw[0] = q[0];
  qw[1] = 0.;
  qw[2] = 0.;
  qw[3] = rGas/gm1*q[0]*bValue[3]; //temperature condition
  qw[4] = 0.;
  rl    = q[0];
  mu    = qa[4];
  p     = .25*Pinv0*max(gamma*mu/(rl*Prn),5.*mu/(3.*rl));
  for (int k=0; k<nq; k++) rb[k] += p*(qw[k]-q[k]);

  //if (viscous) for (int k=0; k<nq; k++) rb[k] -= gv[k];
}
// [BCPenalty]


// [surfaceForces]
void Strand2dFCSPTurbSABcViscousWall::surfaceForces(const double* xs,
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
  mu       = qa [5];
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
