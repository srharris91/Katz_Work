#include "Tri2dFCSPLamBcOutflow.h"


// [Tri2dFCSPLamBcOutflow]
Tri2dFCSPLamBcOutflow::Tri2dFCSPLamBcOutflow()
{
}
// [Tri2dFCSPLamBcOutflow]


// [~Tri2dFCSPLamBcOutflow]
Tri2dFCSPLamBcOutflow::~Tri2dFCSPLamBcOutflow()
{
}
// [~Tri2dFCSPLamBcOutflow]


// [BCVector]
void Tri2dFCSPLamBcOutflow::BCVector(const double* nx,
				     const double* q,
				     const double* qa,
				     double* rb)
{
  double u,v,un,c;
  u     = qa[1];
  v     = qa[2];
  un    = nx[0]*u+nx[1]*v;
  c     = sqrt(gamma*qa[0]/q[0]);
  if (fabs(un/c) < 1.){ //subsonic
    rb[0] = 0.;
    rb[1] = 0.;
    rb[2] = 0.;
    rb[3] = qa[0]-bValue[0]; //pressure
  }
  else{ //supersonic
    rb[0] = 0.;
    rb[1] = 0.;
    rb[2] = 0.;
    rb[3] = 0.;
  }
}
// [BCVector]


// [BCVectorSelfJacobian]
void Tri2dFCSPLamBcOutflow::BCVectorSelfJacobian(const double* nx,
						 const double* q,
						 const double* qa,
						 double* M)
{
  double u,v,qq,un,c;
  u     = qa[1];
  v     = qa[2];
  qq    = u*u+v*v;
  un    = nx[0]*u+nx[1]*v;
  c     = sqrt(gamma*qa[0]/q[0]);

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

  if (fabs(un/c) < 1.){ //subsonic
    M[12] = gm1*qq*.5;
    M[13] =-gm1*u;
    M[14] =-gm1*v;
    M[15] = gm1;
  }
  else{ //supersonic
    M[12] = 0.;
    M[13] = 0.;
    M[14] = 0.;
    M[15] = 0.;
  }
}
// [BCVectorSelfJacobian]


// [BCSelectionMatrix]
void Tri2dFCSPLamBcOutflow::BCSelectionMatrix(const double* nx,
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

  if (fabs(un/c) < 1.){ //subsonic
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
  }
  else{ //supersonic
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
  }
}
// [BCSelectionMatrix]
