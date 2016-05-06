#include "Strand2dFCSPTurbSA.h"


void Strand2dFCSPTurbSA::prepSetup(const int& iPrint0,
				   const int& iTest0,
				   const int& iDebug0,
				   const int& tmp,
				   const int& iSolnFile0,
				   const int& iResdFile0,
				   const int& iErrFile0,
				   int& nq0,
				   int& nqa0,
				   int& ndim0,
				   int& inviscid0,
				   int& viscous0,
				   int& source0,
				   int& sourceMMS0,
				   int& dissipation0,
				   int& nBpatches0,
				   int* iqgrad,
				   int* iqagrad,
				   double* dlim,
				   double* rmsNorm,
				   int& nOutputVars,
				   int* outputVarLength,
				   string* outputVars)
{
  iPrint       = iPrint0;
  iTest        = iTest0;
  iDebug       = iDebug0;
  iSolnFile    = iSolnFile0;
  iResdFile    = iResdFile0;
  iErrFile     = iErrFile0;
  nq0          = nq;
  nqa0         = nqa;
  ndim0        = ndim;
  inviscid0    = inviscid;
  viscous0     = viscous;
  source0      = source;
  sourceMMS0   = sourceMMS;
  dissipation0 = dissipation;
  nBpatches0   = nBpatches;
  nqGradQ      = 5;
  iqgrad [0]   = 1;
  iqgrad [1]   = 1;
  iqgrad [2]   = 1;
  iqgrad [3]   = 1;
  iqgrad [4]   = 1;
  nqaGradQa    = 4;
  iqagrad[0]   = 0; // p
  iqagrad[1]   = 1; // u
  iqagrad[2]   = 1; // v
  iqagrad[3]   = 1; // t
  iqagrad[4]   = 1; // nu
  iqagrad[5]   = 0;
  iqagrad[6]   = 0; 
  iqagrad[7]   = 0;

  int npts = 1;
  int cmp  = 0;
  double r0,p0,u0,v0,t0,e0,q0,qq,nu0;
  double* rValue;
  rValue     = solution.getRefValues();
  p0         = rValue[0];
  u0         = rValue[1];
  v0         = rValue[2];
  t0         = rValue[3];
  nu0        = rValue[4];
  r0         = p0/(rGas*t0);
  qq         = u0*u0+v0*v0;
  e0         = p0/(r0*gm1)+.5*qq;
  q0         = sqrt(qq);
  dlim[0]    = .0625*r0;
  dlim[1]    = .0625*r0*q0;
  dlim[2]    = .0625*r0*q0;
  dlim[3]    = .0625*r0*e0;
  dlim[4]    = .0625*r0*nu0;
  rmsNorm[0] = r0;
  rmsNorm[1] = r0*q0;
  rmsNorm[2] = r0*q0;
  rmsNorm[3] = r0*e0;
  rmsNorm[4] = r0*nu0;

  // SA Constants
  sigma = 2./3.;
  kappa = .41;
  cb1 = .1355;
  cb2 = .622;
  cw1 = cb1/kappa/kappa+(1.+cb2)/sigma;
  cw2 = .3;
  cw3 = 2.;
  cv1 = 7.1;
  cv2 = .7;
  cv3 = .9;
  cn1 = 16.;
  ct3 = 1.2;
  ct4 = 0.5;

  int k=0;
  if (iSolnFile != 0){
    outputVarLength[k] = 1;
    outputVars[k] = "P";
    k++;
    outputVarLength[k] = 1;
    outputVars[k] = "T";
    k++;
    outputVarLength[k] = 1;
    outputVars[k] = "rho";
    k++;
    outputVarLength[k] = 1;
    outputVars[k] = "entropy";
    k++;
    outputVarLength[k] = 1;
    outputVars[k] = "nu_tilde";
    k++;
    outputVarLength[k] = 1;
    outputVars[k] = "eddy_visc";
    k++;
    outputVarLength[k] = 3;
    outputVars[k] = "U";
    k++;
  }

  if (iResdFile != 0){
    outputVarLength[k] = 1;
    outputVars[k] = "R1";
    k++;
    outputVarLength[k] = 1;
    outputVars[k] = "R2";
    k++;
    outputVarLength[k] = 1;
    outputVars[k] = "R3";
    k++;
    outputVarLength[k] = 1;
    outputVars[k] = "R4";
    k++;
    outputVarLength[k] = 1;
    outputVars[k] = "R5";
    k++;
  }

  if (iErrFile != 0){
    outputVarLength[k] = 1;
    outputVars[k] = "E1";
    k++;
    outputVarLength[k] = 1;
    outputVars[k] = "E2";
    k++;
    outputVarLength[k] = 1;
    outputVars[k] = "E3";
    k++;
    outputVarLength[k] = 1;
    outputVars[k] = "E4";
    k++;
    outputVarLength[k] = 1;
    outputVars[k] = "E5";
    k++;
  }
  nOutputVars = k;
}
