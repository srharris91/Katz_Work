#include "Strand2dFCSPLamBcOutflow.h"


// [Strand2dFCSPLamBcOutflow]
Strand2dFCSPLamBcOutflow::Strand2dFCSPLamBcOutflow()
{
}
// [Strand2dFCSPLamBcOutflow]


// [~Strand2dFCSPLamBcOutflow]
Strand2dFCSPLamBcOutflow::~Strand2dFCSPLamBcOutflow()
{
}
// [~Strand2dFCSPLamBcOutflow]


// [BCPenalty]
void Strand2dFCSPLamBcOutflow::BCPenalty(const int& inout,
					 const double* A,
					 const double& Pinv0,
					 const double* q,
					 const double* qa,
					 const double* g,
					 const double* uw,
					 double* rb)
{
  double Ax,Ay,ds,Nx,Ny,Tx,Ty,rl,rul,rvl,re,rqq,p,rhl,rr,rur,rvr,rhr,dd,
    u,v,h,qq,cc,ccr,ccr2,c,ut,un,a,b,m,l[nq],R[nq*nq],S[nq*nq],M[nq*nq];

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
  rr      = g[0];
  rur     = g[1];
  rvr     = g[2];
  re      = g[3];
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

  a       = 0.;
  b       = 0.;
  if (inout > 0){
    if (fabs(un) < c) a = ds*(un+c);
    l[0]    = ds*max(un  ,0.);
    l[1]    = ds*max(un  ,0.);
    l[2]    = ds*max(un+c,0.);
    l[3]    = ds*max(un-c,0.)-a;
  }
  else{
    if (fabs(un) < c) b = ds*(un-c);
    l[0]    = ds*min(un  ,0.);
    l[1]    = ds*min(un  ,0.);
    l[2]    = ds*min(un+c,0.)-b;
    l[3]    = ds*min(un-c,0.);
  }

  R[0 ] = l[0];
  R[1 ] = 0.;
  R[2 ] = l[2]-(a-b);
  R[3 ] = l[3]+(a-b);
  
  R[4 ] = l[0]*u;
  R[5 ] = l[1]*Tx;
  R[6 ] = l[2]*(u+Nx*c)-(a-b)*(u-Nx*c);
  R[7 ] = l[3]*(u-Nx*c)+(a-b)*(u+Nx*c);
  
  R[8 ] = l[0]*v;
  R[9 ] = l[1]*Ty;
  R[10] = l[2]*(v+Ny*c)-(a-b)*(v-Ny*c);
  R[11] = l[3]*(v-Ny*c)+(a-b)*(v+Ny*c);
  
  R[12] = l[0]*qq;
  R[13] = l[1]*ut;
  R[14] = l[2]*(h+un*c)-(a-b)*(h-un*c);
  R[15] = l[3]*(h-un*c)+(a-b)*(h+un*c);
  
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
  
  for (int k=0; k<nq; k++) l[k] = g[k]-q[k];
  
  matmul(nq,nq,nq,&R[0],&S[0],&M[0] );
  matmul(nq,nq,1 ,&M[0],&l[0],&rb[0]);
}
// [BCPenalty]


// [BCPenaltyVis]
void Strand2dFCSPLamBcOutflow::BCPenaltyVis(const double& Pinv0,
					    const double* q,
					    const double* qa,
					    const double* gv,
					    const double* uw,
					    double* rb)
{
  for (int k=0; k<nq; k++) rb[k] =-gv[k];
}
// [BCPenaltyVis]


// [BCPenalty]
void Strand2dFCSPLamBcOutflow::BCPenaltyJacobian(const int& inout,
						 const double* A,
						 const double& Pinv0,
						 const double* q,
						 const double* qa,
						 const double* g,
						 const double* uw,
						 double* M)
{
  double Ax,Ay,ds,Nx,Ny,Tx,Ty,rl,rul,rvl,re,rqq,p,rhl,rr,rur,rvr,rhr,dd,
    u,v,h,qq,cc,ccr,ccr2,c,ut,un,a,b,m,l[nq],R[nq*nq],S[nq*nq];

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
  rr      = g[0];
  rur     = g[1];
  rvr     = g[2];
  re      = g[3];
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

  a       = 0.;
  b       = 0.;
  if (inout > 0){
    if (fabs(un) < c) a = ds*(un+c);
    l[0]    = ds*max(un  ,0.);
    l[1]    = ds*max(un  ,0.);
    l[2]    = ds*max(un+c,0.);
    l[3]    = ds*max(un-c,0.)-a;
  }
  else{
    if (fabs(un) < c) b = ds*(un-c);
    l[0]    = ds*min(un  ,0.);
    l[1]    = ds*min(un  ,0.);
    l[2]    = ds*min(un+c,0.)-b;
    l[3]    = ds*min(un-c,0.);
  }

  R[0 ] = l[0];
  R[1 ] = 0.;
  R[2 ] = l[2]-(a-b);
  R[3 ] = l[3]+(a-b);
  
  R[4 ] = l[0]*u;
  R[5 ] = l[1]*Tx;
  R[6 ] = l[2]*(u+Nx*c)-(a-b)*(u-Nx*c);
  R[7 ] = l[3]*(u-Nx*c)+(a-b)*(u+Nx*c);
  
  R[8 ] = l[0]*v;
  R[9 ] = l[1]*Ty;
  R[10] = l[2]*(v+Ny*c)-(a-b)*(v-Ny*c);
  R[11] = l[3]*(v-Ny*c)+(a-b)*(v+Ny*c);
  
  R[12] = l[0]*qq;
  R[13] = l[1]*ut;
  R[14] = l[2]*(h+un*c)-(a-b)*(h-un*c);
  R[15] = l[3]*(h-un*c)+(a-b)*(h+un*c);
  
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
    
  matmul(nq,nq,nq,&R[0],&S[0],&M[0] );

  for (int k=0; k<nq*nq; k++) M[k] =-M[k];
}
// [BCPenaltyJacobian]
