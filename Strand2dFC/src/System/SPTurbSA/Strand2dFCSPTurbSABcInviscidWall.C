#include "Strand2dFCSPTurbSABcInviscidWall.h"


// [Strand2dFCSPTurbSABcInviscidWall]
Strand2dFCSPTurbSABcInviscidWall::Strand2dFCSPTurbSABcInviscidWall()
{
}
// [Strand2dFCSPTurbSABcInviscidWall]


// [~Strand2dFCSPTurbSABcInviscidWall]
Strand2dFCSPTurbSABcInviscidWall::~Strand2dFCSPTurbSABcInviscidWall()
{
}
// [~Strand2dFCSPTurbSABcInviscidWall]


// [BCPenalty]
void Strand2dFCSPTurbSABcInviscidWall::BCPenalty(const int& inout,
					         const double* A,
					         const double& Pinv0,
					         const double* q,
					         const double* qa,
					         const double* g,
					         const double* uw,
					         double* rb)
{
  double Ax,Ay,ds,Nx,Ny,Tx,Ty,rl,rul,rvl,re,rqq,p,rhl,rr,rur,rvr,rhr,dd,
    rnul,rnur,nu,u,v,h,qq,cc,ccr,ccr2,c,ut,un,m,l[nq],R[nq*nq],S[nq*nq],
    M[nq*nq],qw[nq];

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

  // form the wall state by subtracting off normal velocity
  ut      = rul*Tx+rvl*Ty;
  qw[0]   = q[0];
  qw[1]   = Ny*ut;
  qw[2]   =-Nx*ut;
  qw[3]   = q[3];
  qw[4]   = q[4];

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
  u       =(rul *rl+rur *rr)*dd;
  v       =(rvl *rl+rvr *rr)*dd;
  h       =(rhl *rl+rhr *rr)*dd;
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
}
// [BCPenalty]


// [BCPenaltyVis]
void Strand2dFCSPTurbSABcInviscidWall::BCPenaltyVis(const double& Pinv0,
						 const double* q,
						 const double* qa,
						 const double* gv,
						 const double* uw,
						 double* rb)
{
  for (int k=0; k<nq; k++) rb[k] =-gv[k];
}
// [BCPenaltyVis]


// [BCPenaltyJacobian]
void Strand2dFCSPTurbSABcInviscidWall::BCPenaltyJacobian(const int& inout,
						      const double* A,
						      const double& Pinv0,
						      const double* q,
						      const double* qa,
						      const double* g,
						      const double* uw,
						      double* M)
{
  double Ax,Ay,ds,Nx,Ny,Tx,Ty,rl,rul,rvl,re,rqq,p,rhl,rr,rur,rvr,rhr,dd,
    nu,rnur,rnul,u,v,h,qq,cc,ccr,ccr2,c,ut,un,m,l[nq],R[nq*nq],S[nq*nq],
    P[nq*nq],qw[nq];

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

  // form the wall state by subtracting off normal velocity
  ut      = rul*Tx+rvl*Ty;
  qw[0]   = q[0];
  qw[1]   = Ny*ut;
  qw[2]   =-Nx*ut;
  qw[3]   = q[3];
  qw[4]   = q[4];

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
  u       =(rul *rl+rur *rr)*dd;
  v       =(rvl *rl+rvr *rr)*dd;
  h       =(rhl *rl+rhr *rr)*dd;
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
    
  matmul(nq,nq,nq,&R[0],&S[0],&P[0]);

  S[0 ] = 0.;
  S[1 ] = 0.;
  S[2 ] = 0.;
  S[3 ] = 0.;
  S[4 ] = 0.;

  S[5 ] = 0.;
  S[6 ] = Ny*Tx-1.;
  S[7 ] = Ny*Ty;
  S[8 ] = 0.;
  S[9 ] = 0.;

  S[10] = 0.;
  S[11] =-Nx*Tx;
  S[12] =-Nx*Ty-1.;
  S[13] = 0.;
  S[14] = 0.;

  S[15] = 0.;
  S[16] = 0.;
  S[17] = 0.;
  S[18] = 0.;
  S[19] = 0.;
   
  S[20] = 0.;
  S[21] = 0.; 
  S[22] = 0.;
  S[23] = 0.;
  S[24] = 0.;

  matmul(nq,nq,nq,&P[0],&S[0],&M[0]);

  for (int k=0; k<nq*nq; k++) M[k] *= 2.; //factor of two for well-posedness
}
// [BCPenaltyJacobian]


// [surfaceForces]
void Strand2dFCSPTurbSABcInviscidWall::surfaceForces(const double* xs,
						     const double* ys,
						     const double* q,
						     const double* qa,
						     const double* qx,
						     const double* qy,
						     const double* qax,
						     const double* qay,
						     double* force)
{
  force[0] = qa[0]*ys[0];
  force[1] =-qa[0]*xs[0];
}
// [surfaceForces]
