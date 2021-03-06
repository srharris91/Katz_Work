#include "StrandSPTurbSA.h"
#include <math.h>


void StrandSPTurbSA::lhsSourceJacobian(const int& npts,
				       const double* v,
				       const double* q,
				       const double* qa,
				       const double* qx,
				       const double* qax,
				       double* M)
{
  int iq,iqa,iqax,iqay,iA,iM;
  double rho,rnu,nu,mu,dw,vx,uy,vort,chi,fv1,fv2,S,a1P,a5P,r,g,fw,a,
    a1D,a5D,a1T,a5T;
  
  for (int n=0; n<npts; n++){
    iM       = nq*nq*n;
    iq       = nq *n;
    iqa      = nqa*n;
    iqax     = nqa*n*2;
    iqay     = nqa*n*2+nqa;
    rho      = q[iq  ];
    rnu      = q[iq+4];
    nu       = qa[iqa+4];
    mu       = qa[iqa+5];
    dw       = qa[iqa+7];
    vx       = qax[iqax+2];
    uy       = qax[iqay+1];

    // production
    vort     = fabs(uy-vx);
    chi      = rnu/mu;
    fv1      = pow(chi,3.)/(pow(chi,3.)+pow(cv1,3.));
    fv2      = 1.-chi/(1.+chi*fv1);
    S        = vort+nu*fv2/(kappa*kappa*dw*dw);
    a1P      =-cb1*fv2*nu*nu/(kappa*kappa*dw*dw);
    if (rnu < 0.) a1P = 0.;
    a5P      = cb1*(vort+2.*nu*fv2/(kappa*kappa*dw*dw));
    if (rnu < 0.) a5P = cb1*(1.-ct3)*vort;
    /*
    a1P = 0.;
    a5P      = cb1*S;
    if (rnu < 0.) a5P = cb1*(1.-ct3)*vort;
    */

    // destruction
    r        = nu/(S*kappa*kappa*dw*dw);
    //if (r > 10.) r = 10.;
    g        = r+cw2*(pow(r,6.)-r);
    fw       = g*pow(((1.+pow(cw3,6.))/(pow(g,6.)+pow(cw3,6.))),(1./6.));
    a        = fw;
    if (rnu < 0.) a = -1.;
    a1D      =-cw1*nu*nu/(dw*dw)*a;
    a5D      = cw1*2.*nu/(dw*dw)*a;

    // other term
    a1T      =-cb2*nu*nu/(sigma*v[n]);
    a5T      = cb2*2.*nu/(sigma*v[n]);

    M[iM   ] = 0.;
    M[iM+1 ] = 0.;
    M[iM+2 ] = 0.;
    M[iM+3 ] = 0.;
    M[iM+4 ] = 0.;

    M[iM+5 ] = 0.;
    M[iM+6 ] = 0.;
    M[iM+7 ] = 0.;
    M[iM+8 ] = 0.;
    M[iM+9 ] = 0.;

    M[iM+10] = 0.;
    M[iM+11] = 0.;
    M[iM+12] = 0.;
    M[iM+13] = 0.;
    M[iM+14] = 0.;

    M[iM+15] = 0.;
    M[iM+16] = 0.;
    M[iM+17] = 0.;
    M[iM+18] = 0.;
    M[iM+19] = 0.;

    M[iM+20] = a1P-a1D+a1T;
    M[iM+21] = 0.;
    M[iM+22] = 0.;
    M[iM+23] = 0.;
    M[iM+24] = a5P-a5D+a5T;
  }
}
