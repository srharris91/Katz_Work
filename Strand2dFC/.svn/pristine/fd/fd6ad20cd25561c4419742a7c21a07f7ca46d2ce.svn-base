#include "Strand2dFCSPTurbSA.h"
#include <math.h>

void Strand2dFCSPTurbSA::lhsSourceJacobian(const int& npts,
					   const double* q,
					   const double* qa,
					   const double* qax,
					   const double* qay,
					   double* A)
{
  int iq,iqa,iqax,iA;
  double rho,e,rnu,p,u,v,T,nu,mu,dw,px,py,ux,uy,vx,vy,Tx,Ty,nux,nuy,vort,chi,chi3,
    fv1,fv2,ft2,Sbar,Stil,rr,rr6,g,g6,cw36,fw,dnudq[nq],dmudT,f,dTdq[nq],dmudq[nq],
    dchidq[nq],dfv1dchi,dfv1dq[nq],dfv2dq[nq],dSbardq[nq],dStildSbar,dStildq[nq],
    dft2dq[nq],dDdq[nq],drdq[nq],rtest,dgdq[nq],dgdr,dfwdg,dfwdq[nq];

  for (int n=0; n<npts; n++){
    iq      = nq *n;
    iqa     = nqa*n;
    iqax    = nqa*n;
    iA      = nq*nq*n;
    rho     = q  [iq    ];
    e       = q  [iq  +3]/rho;
    rnu     = q  [iq  +4];
    p       = qa [iqa   ];
    u       = qa [iqa +1];
    v       = qa [iqa +2];
    T       = qa [iqa +3];
    nu      = qa [iqa +4];
    mu      = qa [iqa +5];
    dw      = qa [iqa +7];
    px      = qax[iqax  ];
    py      = qay[iqax  ];
    ux      = qax[iqax+1];
    uy      = qay[iqax+1];
    vx      = qax[iqax+2];
    vy      = qay[iqax+2];
    Tx      = qax[iqax+3];
    Ty      = qay[iqax+3];
    nux     = qax[iqax+4];
    nuy     = qay[iqax+4];

    vort    = fabs(uy-vx);
    chi     = rnu/mu;
    chi3    = chi*chi*chi;
    fv1     = chi3/(chi3+cv1*cv1*cv1);
    fv2     = 1.-chi/(1.+chi*fv1);
    ft2     = ct3*exp(-ct4*chi*chi);
    Sbar    = nu*fv2/(kappa*kappa*dw*dw);
    if (dw < 1.e-13) Sbar = 0.;
    Stil    = vort+Sbar;
    if (Sbar < -cv2*vort)
      Stil  = vort+vort*(vort*cv2*cv2+cv3*Sbar)/((cv3-2.*cv2)*vort-Sbar);
    rr      = nu/(Stil*kappa*kappa*dw*dw);
    rtest   = rr;
    rr      = min(rr,10.);
    rr6     = rr*rr*rr*rr*rr*rr;
    g       = rr+cw2*(rr6-rr);
    g6      = g*g*g*g*g*g;
    cw36    = cw3*cw3*cw3*cw3*cw3*cw3;
    fw      = g*pow(((1.+cw36)/(g6+cw36)),(1./6.));

    dnudq[0] =-nu/rho;
    dnudq[1] = 0.;
    dnudq[2] = 0.;
    dnudq[3] = 0.;
    dnudq[4] = 1./rho;

    transport.getDViscosityDT(1,&p,&T,&dmudT);

    f        = rho*gm1/rGas;
    dTdq[0]  =-f*(e-(u*u+v*v));
    dTdq[1]  =-f*u;
    dTdq[2]  =-f*v;
    dTdq[3]  = f;
    dTdq[4]  = 0.;

    for (int k=0; k<nq; k++) dmudq[k] = dmudT*dTdq[k];

    for (int k=0; k<nq; k++) dchidq[k] =-chi*dmudq[k]/mu;
    dchidq[4] += 1./mu;

    dfv1dchi = 3.*chi*chi/(chi3+cv1*cv1*cv1)*(1.-fv1);

    for (int k=0; k<nq; k++) dfv1dq[k] = dfv1dchi*dchidq[k];

    for (int k=0; k<nq; k++)
      dfv2dq[k] =(chi*(chi*dfv1dq[k]+fv1*dchidq[k])/
		 (1.+chi*fv1)-dchidq[k])/(1.+chi*fv1);
    for (int k=0; k<nq; k++) dSbardq[k] =(nu*dfv2dq[k]+fv2*dnudq[k])/
			                 (kappa*kappa*dw*dw);
    dStildSbar = 1.;
    if (Sbar < -cv2*vort)
      dStildSbar = vort/((cv3-2.*cv2)*vort-Sbar)
	*(cv3+(cv2*cv2*vort+cv3*Sbar)/((cv3-2.*cv2)*vort-Sbar));
    for (int k=0; k<nq; k++) dStildq[k] = dStildSbar*dSbardq[k];
    for (int k=0; k<nq; k++) dft2dq[k] =-2.*ct4*chi*ft2*dchidq[k];
    
    // dDdq
    if (rnu < 0.){
      dDdq[0] = cw1/(dw*dw)*nu*nu;
      dDdq[1] = 0.;
      dDdq[2] = 0.;
      dDdq[3] = 0.;
      dDdq[4] =-cw1/(dw*dw)*2.*nu;
    }
    else{
      for (int k=0; k<nq; k++)
	drdq[k] =(dnudq[k]-nu/Stil*dStildq[k])/(Stil*kappa*kappa*dw*dw);
      if (rtest > 10.) for (int k=0; k<nq; k++) drdq[k] = 0.;
      dgdr    = 1.+cw2*(6.*rr*rr*rr*rr*rr-1.);
      for (int k=0; k<nq; k++) dgdq[k] = dgdr*drdq[k];
      dfwdg   = pow(((1.+cw36)/(g6+cw36)),( 1./6.))
	      - pow(((1.+cw36)/(g6+cw36)),(-5./6.))
	      * ((1.+cw36)/pow((g6+cw36),2))*g6;
      for (int k=0; k<nq; k++) dfwdq[k] = dfwdg*dgdq[k];
      for (int k=0; k<nq; k++)
	dDdq[k] = rho*nu*nu/(dw*dw)*(cw1*dfwdq[k]-cb1/(kappa*kappa)*dft2dq[k]);
      dDdq[0] -= nu*nu*(cw1*fw-cb1/(kappa*kappa)*ft2)/(dw*dw);
      dDdq[4] += 2.*nu*(cw1*fw-cb1/(kappa*kappa)*ft2)/(dw*dw);
    }
    if (dw < 1.e-13) for (int k=0; k<nq; k++) dDdq[k] = 0.;
    
    A[iA   ] = 0.;
    A[iA+1 ] = 0.;
    A[iA+2 ] = 0.;
    A[iA+3 ] = 0.;
    A[iA+4 ] = 0.;
    
    A[iA+5 ] = 0.;
    A[iA+6 ] = 0.;
    A[iA+7 ] = 0.;
    A[iA+8 ] = 0.;
    A[iA+9 ] = 0.;
    
    A[iA+10] = 0.;
    A[iA+11] = 0.;
    A[iA+12] = 0.;
    A[iA+13] = 0.;
    A[iA+14] = 0.;
    
    A[iA+15] = 0.;
    A[iA+16] = 0.;
    A[iA+17] = 0.;
    A[iA+18] = 0.;
    A[iA+19] = 0.;
    
    A[iA+20] =-dDdq[0];
    A[iA+21] =-dDdq[1];
    A[iA+22] =-dDdq[2];
    A[iA+23] =-dDdq[3];
    A[iA+24] =-dDdq[4];	
  }
}
