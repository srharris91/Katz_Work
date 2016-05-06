#include "Strand2dFCSPTurbSA.h"


void Strand2dFCSPTurbSA::rhsDisFlux(const int& npts,
			      const double* Ax,
			      const double* Ay,
			      const double* qL,
			      const double* qR,
			      const double* dq,
			      double* f)
{
  // absolute Jacobian matrix (Roe)
  int iq;
  double ds,Nx,Ny,Tx,Ty,rl,rul,rvl,re,rqq,p,rhl,rr,rur,rvr,rhr,dd,u,v,h,qq,
    cc,ccr,ccr2,c,ut,un,l[nq],a,dlim=.0,R[nq*nq],S[nq*nq],M[nq*nq],rnur,rnul,nu;

  for (int n=0; n<npts; n++){
    iq      = n*nq;

    ds      = sqrt(Ax[n]*Ax[n]+Ay[n]*Ay[n]);
    Nx      = Ax[n]/ds;
    Ny      = Ay[n]/ds;
    Tx      = Ny;
    Ty      =-Nx;

    rl      = qL[iq  ];
    rul     = qL[iq+1];
    rvl     = qL[iq+2];
    re      = qL[iq+3];
    rnul    = qL[iq+4];
    rqq     =(rul*rul+rvl*rvl)/rl;
    p       = gm1*(re-.5*rqq);
    rhl     = re+p;
    rr      = qR[iq  ];
    rur     = qR[iq+1];
    rvr     = qR[iq+2];
    re      = qR[iq+3];
    rnur    = qR[iq+4];
    rqq     =(rur*rur+rvr*rvr)/rr;
    p       = gm1*(re-.5*rqq);
    rhr     = re+p;

    rl      = sqrt(rl);
    rr      = sqrt(rr);
    dd      = 1./(rl+rr);
    rl      = 1./rl;
    rr      = 1./rr;
    u       =(rul *rl+rur *rr)*dd;
    v       =(rvl *rl+rvr *rr)*dd;
    h       =(rhl *rl+rhr *rr)*dd;
    nu      =(rnul*rl+rnur*rr)*dd;
    qq      = .5*(u*u+v*v);
    cc      = gm1*(h-qq);
    ccr     = 1./cc;
    ccr2    = .5*ccr;
    c       = sqrt(cc);
    ut      = u*Tx+v*Ty;
    un      = u*Nx+v*Ny;

    l[0]    = ds*fabs(un  );
    l[1]    = ds*fabs(un  );
    l[2]    = ds*fabs(un+c);
    l[3]    = ds*fabs(un-c);
    l[4]    = ds*fabs(un  );
    a       = ds*dlim*c;
    
    for (int k=0; k<nq; k++) if (l[k] < a) l[k] = .5*(a+l[k]*l[k]/a);

    R[0 ] = l[0];
    R[1 ] = 0.;
    R[2 ] = l[2];
    R[3 ] = l[3];
    R[4 ] = 0.;

    R[5 ] = l[0]*u;
    R[6 ] = l[1]*Tx;
    R[7 ] = l[2]*(u+Nx*c);
    R[8 ] = l[3]*(u-Nx*c);
    R[9 ] = 0.;
    
    R[10] = l[0]*v;
    R[11] = l[1]*Ty;
    R[12] = l[2]*(v+Ny*c);
    R[13] = l[3]*(v-Ny*c);
    R[14] = 0.;

    R[15] = l[0]*qq;
    R[16] = l[1]*ut;
    R[17] = l[2]*(h+un*c);
    R[18] = l[3]*(h-un*c);
    R[19] = 0.;

    R[20] = 0.;
    R[21] = 0.;
    R[22] = l[2]*nu;
    R[23] = l[3]*nu;
    R[24] = l[4];

    S[0 ] =-gm1*ccr*qq+1.;
    S[1 ] = gm1*ccr*u;
    S[2 ] = gm1*ccr*v;
    S[3 ] =-gm1*ccr;
    S[4 ] = 0.;

    S[5 ] =-ut;
    S[6 ] = Tx;
    S[7 ] = Ty;
    S[8 ] = 0.;
    S[9 ] = 0.;

    S[10] = ccr2*(gm1*qq-c*un);
    S[11] =-ccr2*(gm1*u -c*Nx);
    S[12] =-ccr2*(gm1*v -c*Ny);
    S[13] = ccr2* gm1;
    S[14] = 0.;

    S[15] = ccr2*(gm1*qq+c*un);
    S[16] =-ccr2*(gm1*u +c*Nx);
    S[17] =-ccr2*(gm1*v +c*Ny);
    S[18] = ccr2* gm1;
    S[19] = 0.;
 
    S[20] =-nu*gm1*ccr*qq;
    S[21] = nu*gm1*ccr*u;
    S[22] = nu*gm1*ccr*v;
    S[23] =-nu*gm1*ccr;
    S[24] = 1.;


    for (int k=0; k<nq; k++) l[k] = dq[iq+k];

    matmul(nq,nq,nq,&R[0],&S[0],&M[0]);
    matmul(nq,nq,1 ,&M[0],&l[0],&f[0]);
  }
}
