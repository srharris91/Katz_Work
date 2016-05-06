#include "Tri2dFCSPLamBcInflow.h"


// [Tri2dFCSPLamBcInflow]
Tri2dFCSPLamBcInflow::Tri2dFCSPLamBcInflow()
{
}
// [Tri2dFCSPLamBcInflow]


// [~Tri2dFCSPLamBcInflow]
Tri2dFCSPLamBcInflow::~Tri2dFCSPLamBcInflow()
{
}
// [~Tri2dFCSPLamBcInflow]


// [BCVector]
void Tri2dFCSPLamBcInflow::BCVector(const double* nx,
				    const double* q,
				    const double* qa,
				    double* rb)
{
  double p,u,v,t,r,cc,h0,s0,ut0,re,h,s,ut;
  p     = bValue[0];
  u     = bValue[1];
  v     = bValue[2];
  t     = bValue[3];
  r     = p/(rGas*t);
  cc    = gamma*p/r;
  h0    = cc/gm1+.5*(u*u+v*v);
  s0    = p/pow(r,gamma);
  ut0   =-nx[1]*u+nx[0]*v;

  r     = q[0];
  re    = q[3];
  p     = qa[0];
  u     = qa[1];
  v     = qa[2];
  h     =(re+p)/r;
  s     = p/pow(r,gamma);
  ut    =-nx[1]*u+nx[0]*v;

  rb[0] = h -h0 ; //total enthalpy
  rb[1] = s -s0 ; //entropy
  rb[2] = 0.    ;
  rb[3] = ut-ut0; // tangential velocity
}
// [BCVector]


// [BCVectorSelfJacobian]
void Tri2dFCSPLamBcInflow::BCVectorSelfJacobian(const double* nx,
						const double* q,
						const double* qa,
						double* M)
{
  double rr,e,u,v,qq,ut,rrg;
  rr    = 1./q[0];
  e     = rr*q[3];
  u     = qa[1];
  v     = qa[2];
  qq    = u*u+v*v;
  ut    =-nx[1]*u+nx[0]*v;
  rrg   = gm1*pow(rr,gamma);

  M[0 ] = rr*(gm1*qq-gamma*e);
  M[1 ] =-rr*gm1*u;
  M[2 ] =-rr*gm1*v;
  M[3 ] = rr*gamma;

  M[4 ] = rrg*(qq*.5*(gamma+1.)-gamma*e);
  M[5 ] =-rrg*u;
  M[6 ] =-rrg*v;
  M[7 ] = rrg;

  M[8 ] = 0.;
  M[9 ] = 0.;
  M[10] = 0.;
  M[11] = 0.;

  M[12] =-rr*ut;
  M[13] =-rr*nx[1];
  M[14] = rr*nx[0];
  M[15] = 0.;
}
// [BCVectorSelfJacobian]


// [BCSelectionMatrix]
void Tri2dFCSPLamBcInflow::BCSelectionMatrix(const double* nx,
					     const double* q,
					     const double* qa,
					     double* L)
{
  double ccr,qq,u,v,un,c;
  u     = qa[1];
  v     = qa[2];
  un    = nx[0]*u+nx[1]*v;
  qq    = u*u+v*v;
  c     = sqrt(gamma*qa[0]/q[0]);
  ccr   = 1./(c*c);

  L[0 ] = 0.;
  L[1 ] = 0.;
  L[2 ] = 0.;
  L[3 ] = 0.;

  L[4 ] = 0.;
  L[5 ] = 0.;
  L[6 ] = 0.;
  L[7 ] = 0.;

  L[8 ] = ccr*( gm1*qq*.5-c*un   )*.5;
  L[9 ] = ccr*(-gm1*u    +c*nx[0])*.5;
  L[10] = ccr*(-gm1*v    +c*nx[1])*.5;
  L[11] = ccr*( gm1              )*.5;

  L[12] = 0.;
  L[13] = 0.;
  L[14] = 0.;
  L[15] = 0.;
}
// [BCSelectionMatrix]
