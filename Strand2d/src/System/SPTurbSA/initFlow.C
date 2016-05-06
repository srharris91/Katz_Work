#include "StrandSPTurbSA.h"


void StrandSPTurbSA::initFlow(const int& npts,
			   const double* x,
			   double* q)
{
  double* qi = new double[nq*npts];
  solution.getQ(npts,x,qi);

  int i;
  int j=1;
  double p,u,v,t,r,e,qq,nu;
  for (int n=0; n<npts; n++){
    i      = n*nq;
    p      = qi[i  ];
    u      = qi[i+1];
    v      = qi[i+2];
    t      = qi[i+3];
    nu     = qi[i+4];
    r      = p/(rGas*t);
    qq     = u*u+v*v;
    e      = p/(r*gm1)+.5*qq;
    q[i  ] = r;
    q[i+1] = r*u;
    q[i+2] = r*v;
    q[i+3] = r*e;
    q[i+4] = r*nu;
  }

  delete [] qi;
  qi = NULL;
}
