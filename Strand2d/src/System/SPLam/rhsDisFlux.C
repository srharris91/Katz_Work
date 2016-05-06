#include "StrandSPLam.h"


void StrandSPLam::rhsDisFlux(const int& npts,
				const double* A,
				const double* xv,
				const double* ql,
				const double* qr,
				double* f)
{
  // absolute Jacobian matrix (Roe)
  int iq,iA;
  double Ax,Ay,ds,Nx,Ny,Tx,Ty,rl,rul,rvl,re,rqq,p,rhl,rr,rur,rvr,rhr,
    dd,u,v,h,qq,cc,ccr,ccr2,c,ut,un,l[nq],a,dlim=.0,R[nq*nq],S[nq*nq],
    M[nq*nq],eps=1.e-14;

  for (int n=0; n<npts; n++){
    iq      = n*nq;
    iA      = n*ndim;

    Ax      = A[iA  ];
    Ay      = A[iA+1];
    ds      = max(sqrt(Ax*Ax+Ay*Ay),eps);
    Nx      = Ax/ds;
    Ny      = Ay/ds;
    Tx      = Ny;
    Ty      =-Nx;

    rl      = ql[iq  ];
    rul     = ql[iq+1];
    rvl     = ql[iq+2];
    re      = ql[iq+3];
    rqq     =(rul*rul+rvl*rvl)/rl;
    p       = gm1*(re-.5*rqq);
    rhl     = re+p;
    rr      = qr[iq  ];
    rur     = qr[iq+1];
    rvr     = qr[iq+2];
    re      = qr[iq+3];
    rqq     =(rur*rur+rvr*rvr)/rr;
    p       = gm1*(re-.5*rqq);
    rhr     = re+p;

    rl      = sqrt(rl);
    rr      = sqrt(rr);
    dd      = 1./(rl+rr);
    rl      = 1./rl;
    rr      = 1./rr;
    u       =(rul*rl+rur*rr)*dd;
    v       =(rvl*rl+rvr*rr)*dd;
    h       =(rhl*rl+rhr*rr)*dd;
    qq      = .5*(u*u+v*v);
    cc      = gm1*(h-qq);
    ccr     = 1./cc;
    ccr2    = .5*ccr;
    c       = sqrt(cc);
    ut      = u*Tx+v*Ty;
    un      = u*Nx+v*Ny-xv[n]/ds;

    l[0]    = ds*fabs(un  );
    l[1]    = ds*fabs(un  );
    l[2]    = ds*fabs(un+c);
    l[3]    = ds*fabs(un-c);
    a       = dlim*c*ds;
    
    for (int k=0; k<nq; k++) if (l[k] < a) l[k] = .5*(a+l[k]*l[k]/a);

    R[0 ] = l[0];
    R[1 ] = 0.;
    R[2 ] = l[2];
    R[3 ] = l[3];

    R[4 ] = l[0]*u;
    R[5 ] = l[1]*Tx;
    R[6 ] = l[2]*(u+Nx*c);
    R[7 ] = l[3]*(u-Nx*c);

    R[8 ] = l[0]*v;
    R[9 ] = l[1]*Ty;
    R[10] = l[2]*(v+Ny*c);
    R[11] = l[3]*(v-Ny*c);

    R[12] = l[0]*qq;
    R[13] = l[1]*ut;
    R[14] = l[2]*(h+un*c);
    R[15] = l[3]*(h-un*c);
    
    S[0 ] =-gm1*ccr*qq+1.;
    S[1 ] = gm1*ccr*u;
    S[2 ] = gm1*ccr*v;
    S[3 ] =-gm1*ccr;

    S[4 ] =-ut;
    S[5 ] = Tx;
    S[6 ] = Ty;
    S[8 ] = 0.;

    S[8 ] = ccr2*(gm1*qq-c*un);
    S[9 ] =-ccr2*(gm1*u -c*Nx);
    S[10] =-ccr2*(gm1*v -c*Ny);
    S[11] = ccr2* gm1;

    S[12] = ccr2*(gm1*qq+c*un);
    S[13] =-ccr2*(gm1*u +c*Nx);
    S[14] =-ccr2*(gm1*v +c*Ny);
    S[15] = ccr2* gm1;
    
    l[0 ] = .5*(qr[iq  ]-ql[iq  ]);
    l[1 ] = .5*(qr[iq+1]-ql[iq+1]);
    l[2 ] = .5*(qr[iq+2]-ql[iq+2]);
    l[3 ] = .5*(qr[iq+3]-ql[iq+3]);

    matmul(nq,nq,nq,&R[0],&S[0],&M[0 ]);
    matmul(nq,nq,1 ,&M[0],&l[0],&f[iq]);
  }
}
