#include "Strand2dFCSPLam.h"


void Strand2dFCSPLam::initVisFlux(const int& npts,
				  const double* x,
				  double* f)
{
  if (isolution == 3){ //mms

    double* qq = new double[nq*npts];
    double* qx = new double[nq*npts];
    double* qy = new double[nq*npts];
    double* p  = new double[npts];
    double* u  = new double[npts];
    double* ux = new double[npts];
    double* uy = new double[npts];
    double* v  = new double[npts];
    double* vx = new double[npts];
    double* vy = new double[npts];
    double* t  = new double[npts];
    double* tx = new double[npts];
    double* ty = new double[npts];
    double* mu = new double[npts];
    double* k  = new double[npts];

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
      ux[n] = qx[i+1];
      vx[n] = qx[i+2];
      tx[n] = qx[i+3];
      uy[n] = qy[i+1];
      vy[n] = qy[i+2];
      ty[n] = qy[i+3];
    }

    transport.getViscosity   (npts,p,t,mu);
    transport.getConductivity(npts,p,t,k );

    double dd,sxx,syy,sxy,syx;
    for (int n=0; n<npts; n++){
      i      = n*nq*2;
      dd     =(ux[n]+vy[n])/3.;
      sxx    = 2.*mu[n]*(ux[n]-dd);
      syy    = 2.*mu[n]*(vy[n]-dd);
      sxy    =    mu[n]*(uy[n]+vx[n]);
      syx    = sxy;
      f[i  ] = 0.;
      f[i+2] = sxx;
      f[i+4] = sxy;
      f[i+6] = u[n]*sxx+v[n]*sxy+k[n]*tx[n];
      f[i+1] = 0.;
      f[i+3] = syx;
      f[i+5] = syy;
      f[i+7] = u[n]*syx+v[n]*syy+k[n]*ty[n];
    }

    delete [] qq;
    delete [] qx;
    delete [] qy;
    delete [] p;
    delete [] u; 
    delete [] ux;
    delete [] uy;
    delete [] v; 
    delete [] vx;
    delete [] vy;
    delete [] t;
    delete [] tx;
    delete [] ty;
    delete [] mu;
    delete [] k;
  }

  else for (int n=0; n<npts*nq*2; n++) f[n] = 0.;
}
