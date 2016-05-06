#include "Strand2dFCSPTurbSA.h"


void Strand2dFCSPTurbSA::initSource(const int& npts,
				    const double* x,
				    double* s)
{
  if (isolution == 3){ //mms
    for (int n=0; n<npts*nq; n++) s[n] = 0.;

    double* qq  = new double[nq*npts];
    double* qx  = new double[nq*npts];
    double* qy  = new double[nq*npts];
    double* p   = new double[npts];
    double* px  = new double[npts];
    double* py  = new double[npts];
    double* u   = new double[npts];
    double* ux  = new double[npts];
    double* uy  = new double[npts];
    double* v   = new double[npts];
    double* vx  = new double[npts];
    double* vy  = new double[npts];
    double* t   = new double[npts];
    double* tx  = new double[npts];
    double* ty  = new double[npts];
    double* r   = new double[npts];
    double* rx  = new double[npts];
    double* ry  = new double[npts];
    double* h   = new double[npts];
    double* hx  = new double[npts];
    double* hy  = new double[npts];
    double* nu  = new double[npts];
    double* nux = new double[npts];
    double* nuy = new double[npts];

    solution.getQ (npts,x,qq);
    solution.getQx(npts,x,qx);
    solution.getQy(npts,x,qy);

    int i;
    for (int n=0; n<npts; n++){
      i      = n*nq;
      p [n]  = qq[i  ];
      u [n]  = qq[i+1];
      v [n]  = qq[i+2];
      t [n]  = qq[i+3];
      nu[n]  = qq[i+4];

      px[n]  = qx[i  ];
      ux[n]  = qx[i+1];
      vx[n]  = qx[i+2];
      tx[n]  = qx[i+3];
      nux[n] = qx[i+4];

      py[n]  = qy[i  ];
      uy[n]  = qy[i+1];
      vy[n]  = qy[i+2];
      ty[n]  = qy[i+3];
      nuy[n] = qy[i+4]; 
    }
    
    delete [] qq;
    delete [] qx;
    delete [] qy;
    
    state.getDensity  (npts,p,t,r);
    state.getDensityX (npts,p,t,px,tx,rx);
    state.getDensityY (npts,p,t,py,ty,ry);
    state.getEnthalpy (npts,p,t,h);
    state.getEnthalpyX(npts,p,t,px,tx,hx);
    state.getEnthalpyY(npts,p,t,py,ty,hy);
    for (int n=0; n<npts; n++){
      h [n] += (.5*(u[n]*u[n]+v[n]*v[n]));
      hx[n] += (u[n]*ux[n]+v[n]*vx[n]);
      hy[n] += (u[n]*uy[n]+v[n]*vy[n]);
    }

    if (inviscid > 0){
      int is;
      double f1x,f2x,f3x,f4x,f5x,g1y,g2y,g3y,g4y,g5y;
      for (int n=0; n<npts; n++){
	f1x     = r[n]*(                  ux[n])+u[n]*      rx[n];
	f2x     = r[n]*(u[n]*ux [n]+u [n]*ux[n])+u[n]*u [n]*rx[n]+px[n];
	f3x     = r[n]*(u[n]*vx [n]+v [n]*ux[n])+u[n]*v [n]*rx[n];
	f4x     = r[n]*(u[n]*hx [n]+h [n]*ux[n])+u[n]*h [n]*rx[n];
        f5x     = r[n]*(u[n]*nux[n]+nu[n]*ux[n])+u[n]*nu[n]*rx[n];
	
        g1y     = r[n]*(                  vy[n])+v[n]*      ry[n];
	g2y     = r[n]*(v[n]*uy [n]+u [n]*vy[n])+v[n]*u [n]*ry[n];
	g3y     = r[n]*(v[n]*vy [n]+v [n]*vy[n])+v[n]*v [n]*ry[n]+py[n];
	g4y     = r[n]*(v[n]*hy [n]+h [n]*vy[n])+v[n]*h [n]*ry[n];
        g5y     = r[n]*(v[n]*nuy[n]+nu[n]*vy[n])+v[n]*nu[n]*ry[n];
	
        is      = nq*n;
	s[is  ] = f1x+g1y;
	s[is+1] = f2x+g2y;
	s[is+2] = f3x+g3y;
	s[is+3] = f4x+g4y;
        s[is+4] = f5x+g5y;
      }
    }

    if (viscous > 0){
      double* qxx  = new double[nq*npts];
      double* qxy  = new double[nq*npts];
      double* qyy  = new double[nq*npts];
      double* uxx  = new double[npts];
      double* uxy  = new double[npts];
      double* uyy  = new double[npts];
      double* vxx  = new double[npts];
      double* vxy  = new double[npts];
      double* vyy  = new double[npts];
      double* txx  = new double[npts];
      double* txy  = new double[npts];
      double* tyy  = new double[npts];
      double* mu   = new double[npts];
      double* mux  = new double[npts];
      double* muy  = new double[npts];
      double* k    = new double[npts];
      double* kx   = new double[npts];
      double* ky   = new double[npts];
      double* nuxx = new double[npts];
      double* nuyy = new double[npts]; 

      solution.getQxx(npts,x,qxx);
      solution.getQxy(npts,x,qxy);
      solution.getQyy(npts,x,qyy);
      for (int n=0; n<npts; n++){
	i       = n*nq;
	uxx[n]  = qxx[i+1];
	vxx[n]  = qxx[i+2];
	txx[n]  = qxx[i+3];
        nuxx[n] = qxx[i+4];    

	uxy[n]  = qxy[i+1];
	vxy[n]  = qxy[i+2];
	txy[n]  = qxy[i+3];

        uyy[n]  = qyy[i+1];
	vyy[n]  = qyy[i+2];
	tyy[n]  = qyy[i+3];
        nuyy[n] = qyy[i+4];
      }

      delete [] qxx;
      delete [] qxy;
      delete [] qyy;

      transport.getViscosity    (npts,p,t,mu);
      transport.getViscosityX   (npts,p,t,px,tx,mux);
      transport.getViscosityY   (npts,p,t,py,ty,muy);
      transport.getConductivity (npts,p,t,k);
      transport.getConductivityX(npts,p,t,px,tx,kx);
      transport.getConductivityY(npts,p,t,py,ty,ky);

      int is;
      double a,ax,ay,sxx,sxy,syy,sxxx,syyy,sxyx,syxy,qqxx,qqyy;
      double fv1x,fv2x,fv3x,fv4x,fv5x,gv1y,gv2y,gv3y,gv4y,gv5y;
      double rnu,chi,chi3,fv1s,mut,kt,fn,me,ke,eta,chix,chiy,fv1chi,fv1sx,fv1sy,
	fnchi,fnx,fny,etax,etay,mutx,muty,mex,mey,kex,key; 

      for (int n=0; n<npts; n++){        
        rnu      = r[n]*nu[n];
        chi      = rnu/mu[n];
	chi3     = chi*chi*chi;
	fv1s     = chi3/(chi3+cv1*cv1*cv1);
	mut      = rnu*fv1s;
        kt       = mut*rGas*ggm1/PrnT;
        fn       = 1.;
        if (rnu < 0.){
	  mut    = 0.;
	  kt     = 0.;
	  fn     =(cn1+chi3)/(cn1-chi3);
        }
        me       = mu[n]+mut;
        ke       = k[n] +kt;
        eta      = mu[n]+rnu*fn;

	chix     =(rx[n]*nu[n]+r[n]*nux[n]-chi*mux[n])/mu[n];
	chiy     =(ry[n]*nu[n]+r[n]*nuy[n]-chi*muy[n])/mu[n];
        fv1chi   = 3.*chi*chi/(chi3+cv1*cv1*cv1)*(1.-fv1s);
        fv1sx    = fv1chi*chix;
        fv1sy    = fv1chi*chiy;
        fnchi    = 0.;
	if (rnu < 0.) fnchi = 3.*chi*chi/(cn1-chi3)*(1.+fn);
        fnx      = fnchi*chix;
        fny      = fnchi*chiy;
        etax     = mux[n]+rx[n]*nu[n]*fn+r[n]*nux[n]*fn+r[n]*nu[n]*fnx;
        etay     = muy[n]+ry[n]*nu[n]*fn+r[n]*nuy[n]*fn+r[n]*nu[n]*fny; 
        mutx     = rnu*fv1sx+r[n]*nux[n]*fv1s+rx[n]*fv1s*nu[n];
        muty     = rnu*fv1sy+r[n]*nuy[n]*fv1s+ry[n]*fv1s*nu[n]; 
        mex      = mux[n]+mutx;
        mey      = muy[n]+muty;
        kex      = kx[n]+rGas*ggm1/PrnT*mutx;
        key      = ky[n]+rGas*ggm1/PrnT*muty; 
  
        a        =(ux [n]+vy [n])/3.;
	ax       =(uxx[n]+vxy[n])/3.;
	ay       =(uxy[n]+vyy[n])/3.;
	
        sxx      = 2.*me*(ux[n]-a);
	sxy      =    me*(uy[n]+vx[n]);
	syy      = 2.*me*(vy[n]-a);

	sxxx     = 2.*(mex*(ux[n]-a)+me*(uxx[n]-ax));
	syyy     = 2.*(mey*(vy[n]-a)+me*(vyy[n]-ay));
	sxyx     = mex*(uy[n]+vx[n])+me*(uxy[n]+vxx[n]);
	syxy     = mey*(uy[n]+vx[n])+me*(uyy[n]+vxy[n]);
	qqxx     =-(kex*tx[n]+ke*txx[n]);
	qqyy     =-(key*ty[n]+ke*tyy[n]);

        // Viscous Terms
	fv1x     = 0.;
	fv2x     = sxxx;
	fv3x     = sxyx;
	fv4x     = ux[n]*sxx+u[n]*sxxx+vx[n]*sxy+v[n]*sxyx-qqxx;
	fv5x     =(eta*nuxx[n]+nux[n]*etax)/sigma; 

        gv1y     = 0.;
	gv2y     = syxy;
	gv3y     = syyy;
	gv4y     = uy[n]*sxy+u[n]*syxy+vy[n]*syy+v[n]*syyy-qqyy;
	gv5y     =(eta*nuyy[n]+nuy[n]*etay)/sigma; 

        is       = nq*n;
        s[is  ] -=(fv1x+gv1y);
	s[is+1] -=(fv2x+gv2y);
	s[is+2] -=(fv3x+gv3y);
	s[is+3] -=(fv4x+gv4y);
        s[is+4] -=(fv5x+gv5y);
      }

      delete [] uxx;
      delete [] uxy;
      delete [] uyy;
      delete [] vxx;
      delete [] vxy;
      delete [] vyy;
      delete [] txx;
      delete [] txy;
      delete [] tyy;
      delete [] mu; 
      delete [] mux;
      delete [] muy;
      delete [] k;  
      delete [] kx; 
      delete [] ky; 
      delete [] nuxx;
      delete [] nuyy;
    }


    if (source > 0) {
      int is;
      double rnu,chi,chi3,fv1s,mut,kt,fn,vort,fv2,ft2,Sbar,Stil,a,P,rr,rr6,g,g6,
	cw36,fw,D,O,dw;
      double* mu = new double[npts];
      transport.getViscosity(npts,p,t,mu);
      for (int n=0; n<npts; n++){
	dw       = x[2*n+1]; //y-distance above flat plate at y=0.
        rnu      = r[n]*nu[n];
        chi      = rnu/mu[n];
	chi3     = chi*chi*chi;
	fv1s     = chi3/(chi3+cv1*cv1*cv1);
	mut      = rnu*fv1s;
        kt       = mut*rGas*ggm1/PrnT;
        fn       = 1.;
        if (rnu < 0.){
	  mut    = 0.;
	  kt     = 0.;
	  fn     =(cn1+chi3)/(cn1-chi3);
        }

	// production
	vort     = fabs(uy[n]-vx[n]);
	fv2      = 1.-chi/(1.+chi*fv1s);
	ft2      = ct3*exp(-ct4*chi*chi);
	Sbar     = nu[n]*fv2/(kappa*kappa*dw*dw);//0.;
	if (dw < 1.e-13) Sbar = 0.;
	Stil     = vort+Sbar;
	if (Sbar < -cv2*vort)
	  Stil   = vort+vort*(vort*cv2*cv2+cv3*Sbar)/((cv3-2.*cv2)*vort-Sbar);
	a        = Stil*(1.-ft2);
	if (rnu < 0.) a =(1.-ct3)*vort;
	P        = cb1*rnu*a;

	// destruction
	rr       = nu[n]/(Stil*kappa*kappa*dw*dw);
	rr       = min(rr,10.);
	rr6      = rr*rr*rr*rr*rr*rr;
	g        = rr+cw2*(rr6-rr);
	g6       = g*g*g*g*g*g;
	cw36     = cw3*cw3*cw3*cw3*cw3*cw3;
	fw       = g*pow(((1.+cw36)/(g6+cw36)),(1./6.));
	D        = r[n]*(cw1*fw-cb1*ft2/(kappa*kappa))*(nu[n]*nu[n]/(dw*dw));
	if (rnu < 0.) D =-cw1*r[n]*nu[n]*nu[n]/(dw*dw);
	if (dw < 1.e-13) D = 0.;

	// other term
        O        = cb2/sigma*r[n]*(nux[n]*nux[n]+nuy[n]*nuy[n]);
	O       -=(mu[n]+rnu*fn)*(nux[n]*rx[n]+nuy[n]*ry[n])/(r[n]*sigma); 
     
        is       = nq*n;
        s[is  ] -= 0.;
        s[is+1] -= 0.;
        s[is+2] -= 0.;
        s[is+3] -= 0.;     
        s[is+4] -=(P-D+O);
      }
      delete [] mu;
    }

    delete [] p; 
    delete [] px;
    delete [] py;
    delete [] u; 
    delete [] ux;
    delete [] uy;
    delete [] v; 
    delete [] vx;
    delete [] vy;
    delete [] t; 
    delete [] tx;
    delete [] ty;
    delete [] r; 
    delete [] rx;
    delete [] ry;
    delete [] h; 
    delete [] hx;
    delete [] hy;   
    delete [] nu;
    delete [] nux;
    delete [] nuy; 
  }

  else for (int n=0; n<nq*npts; n++) s[n] = 0.;
}
