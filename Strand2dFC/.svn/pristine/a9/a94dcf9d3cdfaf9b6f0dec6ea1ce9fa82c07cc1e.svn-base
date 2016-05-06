#include "Strand2dFCSPTurbSA.h"
#include <math.h>

void Strand2dFCSPTurbSA::rhsSource(const int& npts,
			           const double* q,
			           const double* qa,
			           const double* qax,
			           const double* qay,
			           double* s)
{
  int iq,iqa,iqax;
  double rho,rnu,p,u,v,T,nu,mu,dw,px,py,ux,uy,vx,vy,Tx,Ty,nux,nuy,vort,chi,chi3,
    fv1,fv2,ft2,Sbar,Stil,a,P,rr,rr6,g,g6,cw36,fw,D,rx,ry,fn,O;

  for (int n=0; n<npts; n++){
    iq      = nq *n;
    iqa     = nqa*n;
    iqax    = nqa*n;
    rho     = q  [iq    ];
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

     // production
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
    a       = Stil*(1.-ft2);
    if (rnu < 0.) a =(1.-ct3)*vort;
    P       = cb1*rnu*a;
  
    // destruction
    rr      = nu/(Stil*kappa*kappa*dw*dw);
    rr      = min(rr,10.);
    rr6     = rr*rr*rr*rr*rr*rr;
    g       = rr+cw2*(rr6-rr);
    g6      = g*g*g*g*g*g;
    cw36    = cw3*cw3*cw3*cw3*cw3*cw3;
    fw      = g*pow(((1.+cw36)/(g6+cw36)),(1./6.));
    D       = rho*(cw1*fw-cb1*ft2/(kappa*kappa))*(nu*nu/(dw*dw));
    if (rnu < 0.) D =-cw1*rho*nu*nu/(dw*dw);
    if (dw < 1.e-13) D = 0.;

    // other term
    rx      = px/(rGas*T)-p/(rGas*T*T)*Tx;
    ry      = py/(rGas*T)-p/(rGas*T*T)*Ty;
    fn      = 1.;
    if (rnu < 0.) fn = (cn1+chi3)/(cn1-chi3);
    O       = cb2/sigma*rho*(nux*nux+nuy*nuy);
    O      -=(mu+rnu*fn)*(nux*rx+nuy*ry)/(rho*sigma); 
 
    s[iq  ] = 0.;
    s[iq+1] = 0.;
    s[iq+2] = 0.;
    s[iq+3] = 0.;
    s[iq+4] = P-D+O;	
  }
}
