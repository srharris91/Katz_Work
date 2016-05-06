#include "Tri2dFCSPLam.h"


void Tri2dFCSPLam::rhsVisFluxGalerkin(const int& npts,
				      const double* vr,
				      const double* cx,
				      const double* q,
				      const double* qa,
				      double* f,
				      double* g)
{
  // compute gradients in cell of velocity and temperature
  int iq,iqa,iA,iF;
  double third=1./3.,u[3],v[3],T[3],mu[3],kp[3],Ax[3],Ay[3],ux,uy,vx,vy,Tx,Ty,
    ua,va,mua,kpa,umua,vmua,dd;
  for (int n=0; n<npts; n++){
    iq       = nq *n*3;
    iqa      = nqa*n*3;
    iA       = 2  *n*3;
    iF       = nq *n;

    for (int i=0; i<3; i++){
      u[i]  = qa[iqa+i*nqa+1];
      v[i]  = qa[iqa+i*nqa+2];
      T[i]  = qa[iqa+i*nqa+3];
      mu[i] = qa[iqa+i*nqa+4];
      kp[i] = qa[iqa+i*nqa+5];
      Ax[i] = cx[iA+i*2  ];
      Ay[i] = cx[iA+i*2+1];
    }

    ux =-(Ax[0]*u[0]+Ax[1]*u[1]+Ax[2]*u[2])*vr[n];
    uy =-(Ay[0]*u[0]+Ay[1]*u[1]+Ay[2]*u[2])*vr[n];
    vx =-(Ax[0]*v[0]+Ax[1]*v[1]+Ax[2]*v[2])*vr[n];
    vy =-(Ay[0]*v[0]+Ay[1]*v[1]+Ay[2]*v[2])*vr[n];
    Tx =-(Ax[0]*T[0]+Ax[1]*T[1]+Ax[2]*T[2])*vr[n];
    Ty =-(Ay[0]*T[0]+Ay[1]*T[1]+Ay[2]*T[2])*vr[n];

    ua   = third*(u[0] +u[1] +u[2] );
    va   = third*(v[0] +v[1] +v[2] );
    mua  = third*(mu[0]+mu[1]+mu[2]);
    kpa  = third*(kp[0]+kp[1]+kp[2]);
    umua = third*(u[0]*mu[0]+u[1]*mu[1]+u[2]*mu[2]);
    vmua = third*(v[0]*mu[0]+v[1]*mu[1]+v[2]*mu[2]);
    dd   = third*(ux+vy);

    f[iF  ] = 0.;
    f[iF+1] = 2.*mua*(ux-dd);
    f[iF+2] =    mua*(uy+vx);
    f[iF+3] = 2.*umua*(ux-dd)+vmua*(uy+vx)+kpa*Tx;

    g[iF  ] = 0.;
    g[iF+1] = f[iF+2];
    g[iF+2] = 2.*mua*(vy-dd);
    g[iF+3] = 2.*vmua*(vy-dd)+umua*(uy+vx)+kpa*Ty;
  }
}
