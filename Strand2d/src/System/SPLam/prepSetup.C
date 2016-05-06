#include "StrandSPLam.h"


void StrandSPLam::prepSetup(const int& iPrint0,
			    const int& iTest0,
			    const int& iDebug0,
			    const int& tmp,
			    int& nq0,
			    int& nqa0,
			    int& ndim0,
			    int& inviscid0,
			    int& viscous0,
			    int& source0,
			    int& sourceMMS0,
			    int& dissipation0,
			    int& nBpatches0,
			    int* iqgradT,
			    int* iqagradT,
			    double* dlim,
			    double* rmsNorm)
{
  iPrint       = iPrint0;
  iTest        = iTest0;
  iDebug       = iDebug0;

  nq0          = nq;
  nqa0         = nqa;
  ndim0        = ndim;
  inviscid0    = inviscid;
  viscous0     = viscous;
  source0      = source;
  sourceMMS0   = sourceMMS;
  dissipation0 = dissipation;
  nBpatches0   = nBpatches;
  iqgradT [0]  = 1;
  iqgradT [1]  = 1;
  iqgradT [2]  = 1;
  iqgradT [3]  = 1;
  iqagradT[0]  = 1;
  iqagradT[1]  = 1;
  iqagradT[2]  = 1;
  iqagradT[3]  = 1;
  iqagradT[4]  = 0;
  iqagradT[5]  = 0;

  int npts = 1;
  int cmp  = 0;
  double r0,p0,u0,v0,t0,e0,q0,qq;
  double* rValue;
  rValue       = solution.getRefValues();
  p0           = rValue[0];
  u0           = rValue[1];
  v0           = rValue[2];
  t0           = rValue[3];
  r0           = p0/(rGas*t0);
  qq           = u0*u0+v0*v0;
  e0           = p0/(r0*gm1)+.5*qq;
  q0           = sqrt(qq);
  dlim[0]      = .0625*r0;
  dlim[1]      = .0625*r0*q0;
  dlim[2]      = .0625*r0*q0;
  dlim[3]      = .0625*r0*e0;
  rmsNorm[0]   = r0;
  rmsNorm[1]   = r0*q0;
  rmsNorm[2]   = r0*q0;
  rmsNorm[3]   = r0*e0;
}
