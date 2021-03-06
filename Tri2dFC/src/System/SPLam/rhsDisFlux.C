#include "Tri2dFCSPLam.h"


void Tri2dFCSPLam::rhsDisFlux(const int& npts,
			      const double* Lr0,
			      const double* Ax,
			      const double* Ay,
			      const double* qL,
			      const double* qR,
			      double* f)
{
  // CUSP dissipative flux terms
  int iq;
  double rhL,rhR,rL,rR,dd,uA,vA,hA,cA,unA,mnA,rc0,efix=.5,rc,rp,l,a,
    fdp4,fL,fR,fdp,pL,pR,ruL,rvL,re,rqq,ruR,rvR,qq;

  for (int n=0; n<npts; n++){
    iq      = nq*n;

    // left and right states
    rL      = qL[iq  ];
    ruL     = qL[iq+1];
    rvL     = qL[iq+2];
    re      = qL[iq+3];
    rqq     =(ruL*ruL+rvL*rvL)/rL;
    pL      = gm1*(re-.5*rqq);
    rhL     = re+pL;
    rR      = qR[iq  ];
    ruR     = qR[iq+1];
    rvR     = qR[iq+2];
    re      = qR[iq+3];
    rqq     =(ruR*ruR+rvR*rvR)/rR;
    pR      = gm1*(re-.5*rqq);
    rhR     = re+pR;

    // Roe-averaged state
    rL      = sqrt(rL);
    rR      = sqrt(rR);
    dd      = 1./(rL+rR);
    rL      = 1./rL;
    rR      = 1./rR;
    uA      =(ruL*rL+ruR*rR)*dd;
    vA      =(rvL*rL+rvR*rR)*dd;
    hA      =(rhL*rL+rhR*rR)*dd;
    qq      = .5*(uA*uA+vA*vA);
    cA      = sqrt(gm1*(hA-qq)*(Ax[n]*Ax[n]+Ay[n]*Ay[n]));
    unA     = Ax[n]*uA+Ay[n]*vA;

    // CUSP coefficients
    mnA     = unA/cA;
    rc0     = efix*cA;
    rc      = fabs(unA);
    if (rc < rc0) rc = .5*(rc0+rc*rc/rc0);
    rp      = 1.;
    if (unA < 0.) rp = -1.;
    if (fabs(mnA) < 1.){
      l     = fabs(unA)-cA;
      a     =(fabs(unA)+l)/(fabs(unA)-l);
      rp   *= max(a,0.);
    }
    rc     -= rp*unA;
    unA     =(Ax[n]*qL[iq+1]+Ay[n]*qL[iq+2])/qL[iq];
    fdp4    =-unA*pL;
    fL      = rc+rp*unA;
    unA     =(Ax[n]*qR[iq+1]+Ay[n]*qR[iq+2])/qR[iq];
    fdp4   += unA*pR;
    fdp4   *= rp;
    fR      = rc+rp*unA;
    fdp     = rp*(pR-pL);
    f[iq  ] = fR*qR[iq  ]-fL*qL[iq  ];
    f[iq+1] = fR*qR[iq+1]-fL*qL[iq+1]+Ax[n]*fdp;
    f[iq+2] = fR*qR[iq+2]-fL*qL[iq+2]+Ay[n]*fdp;
    f[iq+3] = fR*qR[iq+3]-fL*qL[iq+3]+fdp4;

  /*
    if (viscous ==-1){
      double tL,tR,muL,muR,pi=4.*atan(1.),alpha=.125,Lr=.5/pi;
      alpha=.0025;
      Lr = Lr0[n];
      tL = pL/(qL[iq]*rGas);
      tR = pR/(qR[iq]*rGas);
      transport.getViscosity(1,&pL,&tL,&muL);
      transport.getViscosity(1,&pR,&tR,&muR);
      a = alpha/Lr*1.4/.75*(muR+muL)/(qR[iq]+qL[iq]);
      for (int k=0; k<nq; k++) f[iq+k] += a*(qR[iq+k]-qL[iq+k]);
    }
  */
  }


  /*
  // absolute Jacobian matrix (Roe)
  int iq;
  double ds,Nx,Ny,Tx,Ty,rl,rul,rvl,re,rqq,p,rhl,rr,rur,rvr,rhr,dd,u,v,h,qq,
    cc,ccr,ccr2,c,ut,un,l[nq],a,dlim=.0,R[nq*nq],S[nq*nq],M[nq*nq];

  for (int n=0; n<npts; n++){
    iq    = n*nq;

    ds    = sqrt(Ax[n]*Ax[n]+Ay[n]*Ay[n]);
    Nx    = Ax[n]/ds;
    Ny    = Ay[n]/ds;
    Tx    = Ny;
    Ty    =-Nx;

    rl    = qL[iq  ];
    rul   = qL[iq+1];
    rvl   = qL[iq+2];
    re    = qL[iq+3];
    rqq   =(rul*rul+rvl*rvl)/rl;
    p     = gm1*(re-.5*rqq);
    rhl   = re+p;
    rr    = qR[iq  ];
    rur   = qR[iq+1];
    rvr   = qR[iq+2];
    re    = qR[iq+3];
    rqq   =(rur*rur+rvr*rvr)/rr;
    p     = gm1*(re-.5*rqq);
    rhr   = re+p;

    rl    = sqrt(rl);
    rr    = sqrt(rr);
    dd    = 1./(rl+rr);
    rl    = 1./rl;
    rr    = 1./rr;
    u     =(rul*rl+rur*rr)*dd;
    v     =(rvl*rl+rvr*rr)*dd;
    h     =(rhl*rl+rhr*rr)*dd;
    qq    = .5*(u*u+v*v);
    cc    = gm1*(h-qq);
    ccr   = 1./cc;
    ccr2  = .5*ccr;
    c     = sqrt(cc);
    ut    = u*Tx+v*Ty;
    un    = u*Nx+v*Ny;

    l[0]  = ds*fabs(un  );
    l[1]  = ds*fabs(un  );
    l[2]  = ds*fabs(un+c);
    l[3]  = ds*fabs(un-c);
    a     = ds*dlim*c;
    
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
    S[7 ] = 0.;

    S[8 ] = ccr2*(gm1*qq-c*un);
    S[9 ] =-ccr2*(gm1*u -c*Nx);
    S[10] =-ccr2*(gm1*v -c*Ny);
    S[11] = ccr2* gm1;

    S[12] = ccr2*(gm1*qq+c*un);
    S[13] =-ccr2*(gm1*u +c*Nx);
    S[14] =-ccr2*(gm1*v +c*Ny);
    S[15] = ccr2* gm1;

    for (int k=0; k<nq; k++) l[k] = qR[iq+k]-qL[iq+k];

    matmul(nq,nq,nq,&R[0],&S[0],&M[0 ]);
    matmul(nq,nq,1 ,&M[0],&l[0],&f[iq]);
  }
  */
}
