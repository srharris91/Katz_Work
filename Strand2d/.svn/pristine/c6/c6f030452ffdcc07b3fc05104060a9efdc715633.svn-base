#include "StrandSPLamBcInflow.h"


// [StrandSPLamBcInflow]
StrandSPLamBcInflow::StrandSPLamBcInflow()
{
}
// [StrandSPLamBcInflow]


// [~StrandSPLamBcInflow]
StrandSPLamBcInflow::~StrandSPLamBcInflow()
{
}
// [~StrandSPLamBcInflow]


// [BCVector]
void StrandSPLamBcInflow::BCVector(const double* nx,
				   const double* wx,
				   const double* qe,
				   const double* qae,
				   const double* q,
				   const double* qa,
				   double* r)
{
  int j=1;
  double Nx,Ny,Tx,Ty,pf,uf,vf,tf,rf,utf,unf,hf,cf,mf;
  double rb,utb,unb,qq,pb,hb,un,du,dv,sb,sf;

  Nx     = nx[0];
  Ny     = nx[1];
  Tx     = nx[1];
  Ty     =-nx[0];

  pf     = bValue[0];
  uf     = bValue[1];
  vf     = bValue[2];
  tf     = bValue[3];
  qq     = uf*uf+vf*vf;
  rf     = pf/(rGas*tf);
  utf    = Tx*uf+Ty*vf;
  unf    = Nx*uf+Ny*vf;
  cf     = sqrt(gamma*rGas*tf);
  hf     = cf*cf/gm1+.5*qq;
  sf     = pf/pow(rf,gamma);

  rb     = q[0];
  pb     = qa[0];
  utb    = Tx*qa[1]+Ty*qa[2];
  unb    = Nx*qa[1]+Ny*qa[2];
  hb     =(q[3]+pb)/rb;
  sb     = pb/pow(rb,gamma);

  mf     = sqrt(qq)/cf;
  un     = unf; //supersonic inflow
  if (mf < 1.) un = Nx*qae[1]+Ny*qae[2]; //subsonic inflow

  //du     = Tx*(utb-utf)+Nx*(unb-un);
  //dv     = Ty*(utb-utf)+Ny*(unb-un);

  du     = Ny*(utb-utf)-Ty*(unb-un);
  dv     =-Nx*(utb-utf)+Tx*(unb-un);

  r[0]   = sb-sf;  //entropy
  r[1]   = du;     //tangential velocity
  r[2]   = dv;     //normal velocity
  r[3]   = hb-hf;  //total enthalpy
}
// [BCVector]


// [BCVectorSelfJacobian]
void StrandSPLamBcInflow::BCVectorSelfJacobian(const double* nx,
					       const double* qe,
					       const double* qae,
					       const double* q,
					       const double* qa,
					       double* M)
{
  double rr,e,rrg,u,v,qq;

  rr    = 1./q[0];
  e     = rr*q[3];
  rrg   = gm1*pow(rr,gamma);
  u     = qa[1];
  v     = qa[2];
  qq    = u*u+v*v;

  M[0 ] = rrg*(qq*.5*(gamma+1.)-gamma*e);
  M[1 ] =-rrg*u;
  M[2 ] =-rrg*v;
  M[3 ] = rrg;

  M[4 ] =-rr*u;
  M[5 ] = rr;
  M[6 ] = 0.;
  M[7 ] = 0.;

  M[8 ] =-rr*v;
  M[9 ] = 0.;
  M[10] = rr;
  M[11] = 0.;

  M[12] = rr*(gm1*qq-gamma*e);
  M[13] =-rr*gm1*u;
  M[14] =-rr*gm1*v;
  M[15] = rr*gamma;
}
// [BCVectorSelfJacobian]


// [BCVectorInteriorJacobian]
void StrandSPLamBcInflow::BCVectorInteriorJacobian(const double* nx,
						   const double* qe,
						   const double* qae,
						   const double* q,
						   const double* qa,
						   double* M)
{
  double uf,vf,qq,tf,cf,mf;
  uf      = bValue[1];
  vf      = bValue[2];
  tf      = bValue[3];
  qq      = uf*uf+vf*vf;
  cf      = sqrt(gamma*rGas*tf);
  mf      = sqrt(qq)/cf;
  for (int n=0; n<nq*nq; n++) M[n] = 0.; //supersonic inflow
  if (mf < 1.){ //subsonic inflow
    double Nx,Ny,Tx,Ty,rr,ue,ve,un;
    Nx    = nx[0];
    Ny    = nx[1];
    Tx    = nx[1];
    Ty    =-nx[0];
    rr    = 1./qe[0];
    ue    = qae[1];
    ve    = qae[2];
    un    = Nx*ue+Ny*ve;
    M[4 ] =-rr*Ty*un;
    M[5 ] = rr*Ty*Nx;
    M[6 ] = rr*Ty*Ny;    
    M[8 ] = rr*Tx*un;
    M[9 ] =-rr*Tx*Nx;
    M[10] =-rr*Tx*Ny; 
  }
}
// [BCVectorInteriorJacobian]
