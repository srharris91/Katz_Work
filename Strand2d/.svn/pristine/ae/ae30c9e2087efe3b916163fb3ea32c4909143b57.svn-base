#include "StrandSPTurbSABcInflow.h"


// [StrandSPTurbSABcInflow]
StrandSPTurbSABcInflow::StrandSPTurbSABcInflow()
{
}
// [StrandSPTurbSABcInflow]


// [~StrandSPTurbSABcInflow]
StrandSPTurbSABcInflow::~StrandSPTurbSABcInflow()
{
}
// [~StrandSPTurbSABcInflow]


// [BCVector]
void StrandSPTurbSABcInflow::BCVector(const double* nx,
				      const double* wx,
				      const double* qe,
				      const double* qae,
				      const double* q,
				      const double* qa,
				      double* r)
{
  /*
  int j=1;
  double Nx,Ny,Tx,Ty,pf,uf,vf,tf,rf,utf,unf,hf,cf,mf;
  double rb,utb,unb,qq,pb,hb,nb,un,du,dv,sb,sf,nf;

  Nx     = nx[0];
  Ny     = nx[1];
  Tx     = nx[1];
  Ty     =-nx[0];

  pf     = bValue[0];
  uf     = bValue[1];
  vf     = bValue[2];
  tf     = bValue[3];
  nf     = bValue[4];
  qq     = uf*uf+vf*vf;
  rf     = pf/(rGas*tf);
  utf    = Tx*uf+Ty*vf;
  unf    = Nx*uf+Ny*vf;
  cf     = sqrt(gamma*rGas*tf);
  hf     = cf*cf/gm1+.5*qq;
  sf     = pf/pow(rf,gamma);

  rb     = q[0];
  pb     = qa[0];
  nb     = qa[4];
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
  r[4]   = nb-nf;  //turbulent eddy viscosity
*/

  int j=1;
  double Nx,Ny,Tx,Ty,Cp,ptf,uf,vf,ttf,nf,qq,utf,unf,cf,mf;
  double rb,pb,tb,nb,utb,unb,ptb,ttb,un,du,dv;

  Nx     = nx[0];
  Ny     = nx[1];
  Tx     = nx[1];
  Ty     =-nx[0];

  Cp     = gamma*rGas/gm1;

  ptf    = bValue[0];
  uf     = bValue[1];
  vf     = bValue[2];
  ttf    = bValue[3];
  nf     = bValue[4];
  qq     = uf*uf+vf*vf;
  utf    = Tx*uf+Ty*vf;
  unf    = Nx*uf+Ny*vf;
  cf     = sqrt(gm1*(Cp*ttf-.5*qq));
  mf     = sqrt(qq)/cf;

  rb     = q[0];
  pb     = qa[0];
  tb     = qa[3];
  nb     = qa[4];
  utb    = Tx*qa[1]+Ty*qa[2];
  unb    = Nx*qa[1]+Ny*qa[2];
  qq     = utb*utb+unb*unb;
  ptb    = pb+.5*rb*qq;
  ttb    = tb+.5*qq/Cp;

  un     = unf; //supersonic inflow
  if (mf < 1.) un = Nx*qae[1]+Ny*qae[2]; //subsonic inflow

  du     = Ny*(utb-utf)-Ty*(unb-un);
  dv     =-Nx*(utb-utf)+Tx*(unb-un);

  r[0]   = ptb-ptf; //total pressure
  r[1]   = du;      //tangential velocity
  r[2]   = dv;      //normal velocity
  r[3]   = ttb-ttf; //total temperature
  r[4]   = nb-nf;   //turbulent eddy viscosity
}
// [BCVector]


// [BCVectorSelfJacobian]
void StrandSPTurbSABcInflow::BCVectorSelfJacobian(const double* nx,
						  const double* qe,
						  const double* qae,
						  const double* q,
						  const double* qa,
						  double* M)
{
  /*
  double rr,e,rrg,u,v,qq,nu;

  rr    = 1./q[0];
  e     = rr*q[3];
  rrg   = gm1*pow(rr,gamma);
  u     = qa[1];
  v     = qa[2];
  nu    = qa[4];
  qq    = u*u+v*v;

  M[0 ] = rrg*(qq*.5*(gamma+1.)-gamma*e);
  M[1 ] =-rrg*u;
  M[2 ] =-rrg*v;
  M[3 ] = rrg;
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

  M[15] = rr*(gm1*qq-gamma*e);
  M[16] =-rr*gm1*u;
  M[17] =-rr*gm1*v;
  M[18] = rr*gamma;
  M[19] = 0.;
  
  M[20] =-rr*nu;
  M[21] = 0.;
  M[22] = 0.;
  M[23] = 0.;
  M[24] = rr;
  */

  double Cp,rr,e,rrg,u,v,qq,nu;

  Cp    = gamma*rGas/gm1;
  rr    = 1./q[0];
  e     = rr*q[3];
  rrg   = rr/Cp;
  u     = qa[1];
  v     = qa[2];
  nu    = qa[4];
  qq    = u*u+v*v;

  M[0 ] =-(2.-gamma)*.5*qq;
  M[1 ] = (2.-gamma)*u;
  M[2 ] = (2.-gamma)*v;
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

  M[15] =-rrg*(gamma*e-gm1*qq);
  M[16] =-rrg*gm1*u;
  M[17] =-rrg*gm1*v;
  M[18] = rrg*gamma;
  M[19] = 0.;
  
  M[20] =-rr*nu;
  M[21] = 0.;
  M[22] = 0.;
  M[23] = 0.;
  M[24] = rr;
}
// [BCVectorSelfJacobian]


// [BCVectorInteriorJacobian]
void StrandSPTurbSABcInflow::BCVectorInteriorJacobian(const double* nx,
						   const double* qe,
						   const double* qae,
						   const double* q,
						   const double* qa,
						   double* M)
{
  /*
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
    M[5 ] =-rr*Ty*un;
    M[6 ] = rr*Ty*Nx;
    M[7 ] = rr*Ty*Ny;    
    M[10] = rr*Tx*un;
    M[11] =-rr*Tx*Nx;
    M[12] =-rr*Tx*Ny; 
  }
  */


  double Cp,uf,vf,qq,ttf,cf,mf;
  Cp     = gamma*rGas/gm1;
  uf     = bValue[1];
  vf     = bValue[2];
  ttf    = bValue[3];
  qq     = uf*uf+vf*vf;
  cf     = sqrt(gm1*(Cp*ttf-.5*qq));
  mf     = sqrt(qq)/cf;
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
    M[5 ] =-rr*Ty*un;
    M[6 ] = rr*Ty*Nx;
    M[7 ] = rr*Ty*Ny;    
    M[10] = rr*Tx*un;
    M[11] =-rr*Tx*Nx;
    M[12] =-rr*Tx*Ny; 
  }
}
// [BCVectorInteriorJacobian]
