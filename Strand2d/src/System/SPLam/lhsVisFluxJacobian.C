#include "StrandSPLam.h"


void StrandSPLam::lhsVisFluxJacobian(const int& npts,
				     const double* A,
				     const double* B,
				     const double* qe,
				     const double* qae,
				     double* M)
{
  int iqa,iA,iM;
  double Ax,Ay,Bx,By,Cxx,Cxy,Cyx,Cyy,a,b,u,v,mu,k,qnA,qnB;
  for (int n=0; n<npts; n++){
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
    u    = qae[iqa+1];
    v    = qae[iqa+2];
    mu   = qae[iqa+4];
    k    = qae[iqa+5];
    qnA  = Ax*u+Ay*v;
    qnB  = Bx*u+By*v;


    /*
    mu = 1.634156530983452e-5;
    k  = mu*1.4*287.06/(.75*.4);
    u  = 6.94449998200014e1;
    v  = 0.;
    */


    M[iM   ] = 0.;
    M[iM+1 ] = 0.;
    M[iM+2 ] = 0.;
    M[iM+3 ] = 0.;

    M[iM+4 ] = 0.;
    M[iM+5 ] = mu*( a+Cxx);
    M[iM+6 ] = mu*(-b+Cxy);
    M[iM+7 ] = 0.;

    M[iM+8 ] = 0.;
    M[iM+9 ] = mu*( b+Cyx);
    M[iM+10] = mu*( a+Cyy);
    M[iM+11] = 0.;

    M[iM+12] = 0.;
    M[iM+13] = mu*(a*u+Ax*qnB-2./3.*Bx*qnA);
    M[iM+14] = mu*(a*v+Ay*qnB-2./3.*By*qnA);
    M[iM+15] = a*k;


    /*
    double sr;
    sr = fabs(mu*(2.*a+Cxx+Cyy));
    if (fabs(a*k) > sr) sr = fabs(a*k);

    M[iM   ] = 0.;
    M[iM+1 ] = 0.;
    M[iM+2 ] = 0.;
    M[iM+3 ] = 0.;

    M[iM+4 ] = 0.;
    M[iM+5 ] = sr;
    M[iM+6 ] = 0.;
    M[iM+7 ] = 0.;

    M[iM+8 ] = 0.;
    M[iM+9 ] = 0.;
    M[iM+10] = sr;
    M[iM+11] = 0.;

    M[iM+12] = 0.;
    M[iM+13] = 0.;
    M[iM+14] = 0.;
    M[iM+15] = sr;
    */
  }
}


void StrandSPLam::lhsVisFluxJacobian(const int& npts,
				     const double* A,
				     const double* q,
				     const double* qa,
				     const double* qx,
				     const double* qy,
				     const double* qax,
				     const double* qay,
				     double* M)
{
  int iq,iqa,iA,iM;
  double Ax,Ay,r,rr,rrr,a,u,v,qq,mu,k,rx,ry,rex,rey,ux,uy,vx,vy,dd,sxx,sxy,syy,
    re,e,dudq1,dudq2,dvdq1,dvdq3,
    duxdq1,duxdq2,duydq1,duydq2,dvxdq1,dvxdq3,dvydq1,dvydq3,
    dtxdq1,dtxdq2,dtxdq3,dtxdq4,dtydq1,dtydq2,dtydq3,dtydq4,
    dsxxdq1,dsxxdq2,dsxxdq3,dsxydq1,dsxydq2,dsxydq3,dsyydq1,dsyydq2,dsyydq3,
    dqxdq1,dqxdq2,dqxdq3,dqxdq4,dqydq1,dqydq2,dqydq3,dqydq4;

  for (int n=0; n<npts; n++){
    iq   = nq*n;
    iqa  = nqa*n;
    iA   = ndim*n;
    iM   = nq*nq*n;
    Ax   = A[iA  ];
    Ay   = A[iA+1];
    r    = q  [iq   ];
    re   = q  [iq+3 ];
    rr   = 1./r;
    rrr  = rr/r;
    e    = re*rr;
    a    = gm1*rrr/rGas;
    u    = qa [iqa+1];
    v    = qa [iqa+2];
    qq   = u*u+v*v;
    rx   = qx [iq   ];
    ry   = qy [iq   ];
    rex  = qx [iq+3 ];
    rey  = qy [iq+3 ];
    mu   = qa [iqa+4];
    k    = qa [iqa+5];
    ux   = qax[iqa+1];
    uy   = qay[iqa+1];
    vx   = qax[iqa+2];
    vy   = qay[iqa+2];
    dd   =(ux+vy)/3.;
    sxx  = 2.*mu*(ux-dd);
    sxy  =    mu*(uy+vx);
    syy  = 2.*mu*(vy-dd);

    dudq1 =-u*rr;
    dudq2 = rr;
    dvdq1 =-v*rr;
    dvdq3 = rr;

    duxdq1 = rrr*(u*rx-r*ux);
    duxdq2 =-rrr*rx;

    duydq1 = rrr*(u*ry-r*uy);
    duydq2 =-rrr*ry;

    dvxdq1 = rrr*(v*rx-r*vx);
    dvxdq3 =-rrr*rx;

    dvydq1 = rrr*(v*ry-r*vy);
    dvydq3 =-rrr*ry;

    dtxdq1 = a*(2.*r*(u*ux+v*vx)+rx*(2.*e-qq)-rex);
    dtxdq2 =-a*(r*ux-u*rx);
    dtxdq3 =-a*(r*vx-v*rx);
    dtxdq4 =-a*rx;

    dtydq1 = a*(2.*r*(u*uy+v*vy)+ry*(2.*e-qq)-rey);
    dtydq2 =-a*(r*uy-u*ry);
    dtydq3 =-a*(r*vy-v*ry);
    dtydq4 =-a*ry;

    dsxxdq1 = 2./3.*mu*(2.*duxdq1-dvydq1);
    dsxxdq2 = 2./3.*mu*(2.*duxdq2       );
    dsxxdq3 = 2./3.*mu*(         -dvydq3);

    dsxydq1 = mu*(duydq1+dvxdq1);
    dsxydq2 = mu*(duydq2       );
    dsxydq3 = mu*(       dvxdq3);

    dsyydq1 = 2./3.*mu*(2.*dvydq1-duxdq1);
    dsyydq2 = 2./3.*mu*(         -duxdq2);
    dsyydq3 = 2./3.*mu*(2.*dvydq3       );

    dqxdq1 =-k*dtxdq1;
    dqxdq2 =-k*dtxdq2;
    dqxdq3 =-k*dtxdq3;
    dqxdq4 =-k*dtxdq4;

    dqydq1 =-k*dtydq1;
    dqydq2 =-k*dtydq2;
    dqydq3 =-k*dtydq3;
    dqydq4 =-k*dtydq4;

    M[iM   ] = 0.;
    M[iM+1 ] = 0.;
    M[iM+2 ] = 0.;
    M[iM+3 ] = 0.;

    M[iM+4 ] = Ax*dsxxdq1+Ay*dsxydq1;
    M[iM+5 ] = Ax*dsxxdq2+Ay*dsxydq2;
    M[iM+6 ] = Ax*dsxxdq3+Ay*dsxydq3;
    M[iM+7 ] = 0.;

    M[iM+8 ] = Ax*dsxydq1+Ay*dsyydq1;
    M[iM+9 ] = Ax*dsxydq2+Ay*dsyydq2;
    M[iM+10] = Ax*dsxydq3+Ay*dsyydq3;
    M[iM+11] = 0.;

    M[iM+12] = dudq1*(Ax*sxx+Ay*sxy)+u*M[iM+4 ]
             + dvdq1*(Ax*sxy+Ay*syy)+v*M[iM+8 ]
             - Ax*dqxdq1-Ay*dqydq1;
    M[iM+13] = dudq2*(Ax*sxx+Ay*sxy)+u*M[iM+5 ]
             + v*M[iM+9 ]
             - Ax*dqxdq2-Ay*dqydq2;
    M[iM+14] = u*M[iM+6 ]
             + dvdq3*(Ax*sxy+Ay*syy)+v*M[iM+10]
             - Ax*dqxdq3-Ay*dqydq3;
    M[iM+15] =-Ax*dqxdq4-Ay*dqydq4;
  }
}
