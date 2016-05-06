#include "Strand2dFCSPTurbSA.h"


void Strand2dFCSPTurbSA::initVisFlux(const int& npts,
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
    double* nu = new double[npts];
    double* nux= new double[npts];
    double* nuy= new double[npts];

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
      nu[n] = qq[i+4];
      
      ux[n] = qx[i+1];
      vx[n] = qx[i+2];
      tx[n] = qx[i+3];
      nux[n]= qx[i+4];
      
      uy[n] = qy[i+1];
      vy[n] = qy[i+2];
      ty[n] = qy[i+3];
      nuy[n]= qy[i+4];
    }

    transport.getViscosity   (npts,p,t,mu);
    transport.getConductivity(npts,p,t,k );

    double rho,rnu,chi,chi3,fv1,mut,kt,me,ke,fn,dd,sxx,syy,sxy,syx,qxx,qyy,nvx,nvy;
    for (int n=0; n<npts; n++){
      i      = n*nq*2;
      rho    = p[n]/(rGas*t[n]);
      rnu    = rho*nu[n];
      chi    = rnu/mu[n];
      chi3   = chi*chi*chi;
      fv1    = chi3/(chi3+cv1*cv1*cv1);
      mut    = rnu*fv1;
      kt     = mut*rGas*ggm1/PrnT;
      fn     =(mu[n]+rnu)/sigma;
      if (rnu < 0.){
	fn   =(cn1+chi3)/(cn1-chi3);
	fn   =(mu[n]+rnu*fn)/sigma;
	mut  = 0.;
	kt   = 0.;
      }
      me     = mu[n]+mut;
      ke     = k[n] +kt;

      dd     =(ux[n]+vy[n])/3.;
      sxx    = 2.*me*(ux[n]-dd);
      syy    = 2.*me*(vy[n]-dd);
      sxy    =    me*(uy[n]+vx[n]);
      syx    = sxy;
      qxx    =-ke*tx[n];
      qyy    =-ke*ty[n];
      nvx    = fn*nux[n];
      nvy    = fn*nuy[n];

      f[i  ] = 0.;
      f[i+2] = sxx;
      f[i+4] = sxy;
      f[i+6] = u[n]*sxx+v[n]*sxy-qxx;
      f[i+8] = nvx;
     
      f[i+1] = 0.;
      f[i+3] = syx;
      f[i+5] = syy;
      f[i+7] = u[n]*syx+v[n]*syy-qyy;
      f[i+9] = nvy;

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
    delete [] nu;
    delete [] nux;
    delete [] nuy;
  }

  else for (int n=0; n<npts*nq*2; n++) f[n] = 0.;


}
