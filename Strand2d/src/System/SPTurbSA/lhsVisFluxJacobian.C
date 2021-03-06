#include "StrandSPTurbSA.h"


void StrandSPTurbSA::lhsVisFluxJacobian(const int& npts,
					const double* A,
					const double* B,
					const double* qe,
					const double* qae,
					double* M)
{
  int iq,iqa,iA,iM;
  double Ax,Ay,Bx,By,Cxx,Cxy,Cyx,Cyy,a,b,rnu,u,v,mu,k,nu,qnA,qnB,
    fn,chi,fv1,mut,kt;
  
  for (int n=0; n<npts; n++){
    iq   = nq*n;
    iqa  = nqa*n;
    iA   = ndim*n;
    iM   = nq*nq*n;
    Ax   = A[iA  ];
    Ay   = A[iA+1];
    Bx   = B[iA  ];
    By   = B[iA+1];
    Cxx  = Ax*Bx/3.;
    Cxy  = Ax*By/3.;
    Cyx  = Ay*Bx/3.;
    Cyy  = Ay*By/3.;
    a    = Ax*Bx+Ay*By;
    b    = Ax*By-Ay*Bx;
    rnu  = qe[iq+4];
    u    = qae[iqa+1];
    v    = qae[iqa+2];
    nu   = qae[iqa+4];
    mu   = qae[iqa+5];
    k    = qae[iqa+6];
    qnA  = Ax*u+Ay*v;
    qnB  = Bx*u+By*v;

    chi  = rnu/mu;
    fv1  = pow(chi,3.)/(pow(chi,3.)+pow(cv1,3.));
    mut  = fv1*rnu;
    kt   = mut*rGas*ggm1/prnT;
    fn   =(mu+rnu)/sigma;
    if (rnu < 0.){
      fn  =(cn1+pow(chi,3.))/(cn1-pow(chi,3.));
      fn  =(mu+rnu*fn)/sigma;
      mut = 0.;
      kt  = 0.;
    }
    mu  += mut;
    k   += kt;

    M[iM   ] = 0.;
    M[iM+1 ] = 0.;
    M[iM+2 ] = 0.;
    M[iM+3 ] = 0.;
    M[iM+4 ] = 0.;

    M[iM+5 ] = 0.;
    M[iM+6 ] = mu*( a+Cxx);
    M[iM+7 ] = mu*(-b+Cxy);
    M[iM+8 ] = 0.;
    M[iM+9 ] = 0.;

    M[iM+10] = 0.;
    M[iM+11] = mu*( b+Cyx);
    M[iM+12] = mu*( a+Cyy);
    M[iM+13] = 0.;
    M[iM+14] = 0.;

    M[iM+15] = 0.;
    M[iM+16] = mu*(a*u+Ax*qnB-2./3.*Bx*qnA);
    M[iM+17] = mu*(a*v+Ay*qnB-2./3.*By*qnA);
    M[iM+18] = a*k;
    M[iM+19] = 0.;
    
    M[iM+20] = 0.;
    M[iM+21] = 0.;
    M[iM+22] = 0.;
    M[iM+23] = 0.;
    M[iM+24] = fn*a;
  }
}
