#include "StrandSPTurbSA.h"
#include <math.h>

void StrandSPTurbSA::rhsSource(const int& npts,
			       const double* q,
			       const double* qa,
			       const double* qx,
			       const double* qax,
			       double* r)
{
  int iq,iqa,iqax,iqay;
  double rho,rnu,nu,mu,dw,vx,uy,nux,nuy,vort,chi,fv1,fv2,S,a,rr,g,fw,P,D,T;

  for (int n=0; n<npts; n++){
    iq      = nq *n;
    iqa     = nqa*n;
    iqax    = nqa*n*2;
    iqay    = nqa*n*2+nqa;
    rho     = q[iq  ];
    rnu     = q[iq+4];
    nu      = qa[iqa+4];
    mu      = qa[iqa+5];
    dw      = qa[iqa+7];
    vx      = qax[iqax+2];
    uy      = qax[iqay+1];
    nux     = qax[iqax+4];
    nuy     = qax[iqay+4];

    // production
    vort    = fabs(uy-vx);
    chi     = rnu/mu;
    fv1     = pow(chi,3.)/(pow(chi,3.)+pow(cv1,3.));
    fv2     = 1.-chi/(1.+chi*fv1);
    S       = vort+nu*fv2/(kappa*kappa*dw*dw);
    a       = S;
    if (rnu < 0.) a =(1.-ct3)*vort;
    P       = cb1*rnu*a;

    // destruction
    rr      = nu/(S*kappa*kappa*dw*dw);
    //if (rr > 10.) rr = 10.;
    g       = rr+cw2*(pow(rr,6.)-rr);
    fw      = g*pow(((1.+pow(cw3,6.))/(pow(g,6.)+pow(cw3,6.))),(1./6.));
    a       = fw;
    if (rnu < 0.) a = -1.;
    D       = cw1*rho*nu*nu/(dw*dw)*a;

    // other term
    T       = cb2/sigma*rho*(nux*nux+nuy*nuy);

    r[iq  ] = 0.;
    r[iq+1] = 0.;
    r[iq+2] = 0.;
    r[iq+3] = 0.;
    r[iq+4] = P-D+T;
  }
}
