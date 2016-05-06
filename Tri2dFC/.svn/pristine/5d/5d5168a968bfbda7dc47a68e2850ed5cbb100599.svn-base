#include "Tri2dFCSPLam.h"


void Tri2dFCSPLam::initSource(const int& npts,
			      const double* x,
			      double* s)
{
  if (isolution == 3){ //mms
    for (int n=0; n<npts*nq; n++) s[n] = 0.;

    double* qq = new double[nq*npts];
    double* qx = new double[nq*npts];
    double* qy = new double[nq*npts];
    double* p  = new double[npts];
    double* px = new double[npts];
    double* py = new double[npts];
    double* u  = new double[npts];
    double* ux = new double[npts];
    double* uy = new double[npts];
    double* v  = new double[npts];
    double* vx = new double[npts];
    double* vy = new double[npts];
    double* t  = new double[npts];
    double* tx = new double[npts];
    double* ty = new double[npts];
    double* r  = new double[npts];
    double* rx = new double[npts];
    double* ry = new double[npts];
    double* h  = new double[npts];
    double* hx = new double[npts];
    double* hy = new double[npts];

    solution.getQ (npts,x,qq);
    solution.getQx(npts,x,qx);
    solution.getQy(npts,x,qy);

    int i;
    for (int n=0; n<npts; n++){
      i     = n*nq;
      p [n] = qq[i  ];
      u [n] = qq[i+1];
      v [n] = qq[i+2];
      t [n] = qq[i+3];
      px[n] = qx[i  ];
      ux[n] = qx[i+1];
      vx[n] = qx[i+2];
      tx[n] = qx[i+3];
      py[n] = qy[i  ];
      uy[n] = qy[i+1];
      vy[n] = qy[i+2];
      ty[n] = qy[i+3];
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
      double f1x,f2x,f3x,f4x,g1y,g2y,g3y,g4y;
      for (int n=0; n<npts; n++){
	f1x     = r[n]*(                ux[n])+u[n]*     rx[n];
	f2x     = r[n]*(u[n]*ux[n]+u[n]*ux[n])+u[n]*u[n]*rx[n]+px[n];
	f3x     = r[n]*(u[n]*vx[n]+v[n]*ux[n])+u[n]*v[n]*rx[n];
	f4x     = r[n]*(u[n]*hx[n]+h[n]*ux[n])+u[n]*h[n]*rx[n];
	g1y     = r[n]*(                vy[n])+v[n]*     ry[n];
	g2y     = r[n]*(v[n]*uy[n]+u[n]*vy[n])+v[n]*u[n]*ry[n];
	g3y     = r[n]*(v[n]*vy[n]+v[n]*vy[n])+v[n]*v[n]*ry[n]+py[n];
	g4y     = r[n]*(v[n]*hy[n]+h[n]*vy[n])+v[n]*h[n]*ry[n];
	is      = nq*n;
	s[is  ] = f1x+g1y;
	s[is+1] = f2x+g2y;
	s[is+2] = f3x+g3y;
	s[is+3] = f4x+g4y;
      }
    }

    if (viscous > 0){
      double* qxx = new double[nq*npts];
      double* qxy = new double[nq*npts];
      double* qyy = new double[nq*npts];
      double* uxx = new double[npts];
      double* uxy = new double[npts];
      double* uyy = new double[npts];
      double* vxx = new double[npts];
      double* vxy = new double[npts];
      double* vyy = new double[npts];
      double* txx = new double[npts];
      double* txy = new double[npts];
      double* tyy = new double[npts];
      double* mu  = new double[npts];
      double* mux = new double[npts];
      double* muy = new double[npts];
      double* k   = new double[npts];
      double* kx  = new double[npts];
      double* ky  = new double[npts];

      solution.getQxx(npts,x,qxx);
      solution.getQxy(npts,x,qxy);
      solution.getQyy(npts,x,qyy);
      for (int n=0; n<npts; n++){
	i      = n*nq;
	uxx[n] = qxx[i+1];
	vxx[n] = qxx[i+2];
	txx[n] = qxx[i+3];
	uxy[n] = qxy[i+1];
	vxy[n] = qxy[i+2];
	txy[n] = qxy[i+3];
	uyy[n] = qyy[i+1];
	vyy[n] = qyy[i+2];
	tyy[n] = qyy[i+3];
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
      double fv1x,fv2x,fv3x,fv4x,gv1y,gv2y,gv3y,gv4y;
      for (int n=0; n<npts; n++){
	a        =(ux [n]+vy [n])/3.;
	ax       =(uxx[n]+vxy[n])/3.;
	ay       =(uxy[n]+vyy[n])/3.;
	sxx      = 2.*mu[n]*(ux[n]-a);
	sxy      =    mu[n]*(uy[n]+vx[n]);
	syy      = 2.*mu[n]*(vy[n]-a);
	sxxx     = 2.*(mux[n]*(ux[n]-a)+mu[n]*(uxx[n]-ax));
	syyy     = 2.*(muy[n]*(vy[n]-a)+mu[n]*(vyy[n]-ay));
	sxyx     = mux[n]*(uy[n]+vx[n])+mu[n]*(uxy[n]+vxx[n]);
	syxy     = muy[n]*(uy[n]+vx[n])+mu[n]*(uyy[n]+vxy[n]);
	qqxx     =-(kx[n]*tx[n]+k[n]*txx[n]);
	qqyy     =-(ky[n]*ty[n]+k[n]*tyy[n]);
	fv1x     = 0.;
	fv2x     = sxxx;
	fv3x     = sxyx;
	fv4x     = ux[n]*sxx+u[n]*sxxx+vx[n]*sxy+v[n]*sxyx-qqxx;
	gv1y     = 0;
	gv2y     = syxy;
	gv3y     = syyy;
	gv4y     = uy[n]*sxy+u[n]*syxy+vy[n]*syy+v[n]*syyy-qqyy;
	is       = nq*n;
	s[is  ] -=(fv1x+gv1y);
	s[is+1] -=(fv2x+gv2y);
	s[is+2] -=(fv3x+gv3y);
	s[is+3] -=(fv4x+gv4y);
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
  }

  else for (int n=0; n<nq*npts; n++) s[n] = 0.;
}
