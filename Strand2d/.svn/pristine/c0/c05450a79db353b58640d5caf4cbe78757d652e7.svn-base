#include "StrandSPTurbSABcOutflow.h"


// [StrandSPTurbSABcOutflow]
StrandSPTurbSABcOutflow::StrandSPTurbSABcOutflow()
{
}
// [StrandSPTurbSABcOutflow]


// [~StrandSPTurbSABcOutflow]
StrandSPTurbSABcOutflow::~StrandSPTurbSABcOutflow()
{
}
// [~StrandSPTurbSABcOutflow]


// [BCVector]
void StrandSPTurbSABcOutflow::BCVector(const double* nx,
				       const double* wx,
				       const double* qe,
				       const double* qae,
				       const double* q,
				       const double* qa,
				       double* r)
{
  int j=1;
  double pe,ue,ve,te,nue,ce,me,p;

  pe     = qae[0];
  ue     = qae[1];
  ve     = qae[2];
  te     = qae[3];
  nue    = qae[4];
  ce     = gamma*rGas*te;
  me     = sqrt((ue*ue+ve*ve)/ce);

  p      = bValue[0]; // subsonic outflow
  if (me > 1.) p = pe; //supersonic outflow

  r[0] = qa[0]-p;
  r[1] = qa[1]-ue;
  r[2] = qa[2]-ve;
  r[3] = qa[3]-te;
  r[4] = qa[4]-nue;
}
// [BCVector]


// [BCVectorSelfJacobian]
void StrandSPTurbSABcOutflow::BCVectorSelfJacobian(const double* nx,
						   const double* qe,
						   const double* qae,
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

  M[0 ] = gm1*qq*.5;
  M[1 ] =-gm1*u;
  M[2 ] =-gm1*v;
  M[3 ] = gm1;
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


// [BCVectorInteriorJacobian]
void StrandSPTurbSABcOutflow::BCVectorInteriorJacobian(const double* nx,
						       const double* qe,
						       const double* qae,
						       const double* q,
						       const double* qa,
						       double* M)
{
  double rr,e,u,v,t,c,m,qq,grr,nu;
  rr    = 1./qe[0];
  e     = qe[3]*rr;
  u     = qae[1];
  v     = qae[2];
  t     = qae[3];
  nu    = qae[4];
  c     = gamma*rGas*t;
  m     = sqrt((u*u+v*v)/c);
  qq    = u*u+v*v;
  grr   = gm1*rr/rGas;

  if (m > 1.){ //supersonic outflow
    M[0 ] =-gm1*qq*.5;
    M[1 ] = gm1*u;
    M[2 ] = gm1*v;
    M[3 ] =-gm1;
  }
  else{ // subsonic outflow
    M[0 ] = 0.;
    M[1 ] = 0.;
    M[2 ] = 0.;
    M[3 ] = 0.;
  }
  
  M[4 ] = 0.;
 
  M[5 ] = rr*u;
  M[6 ] =-rr;
  M[7 ] = 0.;
  M[8 ] = 0.;
  M[9 ] = 0.;

  M[10] = rr*v;
  M[11] = 0.;
  M[12] =-rr;
  M[13] = 0.;
  M[14] = 0.;

  M[15] =-grr*(qq-e);
  M[16] = grr*u;
  M[17] = grr*v;
  M[18] =-grr;
  M[19] = 0.;
  
  M[20] = rr*nu;
  M[21] = 0.;
  M[22] = 0.;
  M[23] = 0.;
  M[24] =-rr;
}
// [BCVectorInteriorJacobian]
