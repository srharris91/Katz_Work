#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::hessian(const int& nqp,
				 const double* p,
				 double* pxx)
{
  // FEM Hessian
  int k1x,k1y,k1xx,k1xy,k1yy,k2,k3;
  double dx,dy,dxx,dxy,dyy;
  for (int n=0; n<nNode*nqp*5; n++) pxx[n] = 0.;
  for (int n=0; n<nNode; n++){
    k1x  = n*nqp*5;
    k1y  = k1x+nqp;
    k1xx = k1y+nqp;
    k1xy = k1xx+nqp;
    k1yy = k1xy+nqp;
    for(int i=psp2(n); i<psp2(n+1); i++){
      k2  = psp1(i)*nqp;
      dx  = gx (i,0);
      dy  = gx (i,1);
      dxx = gxx(i,0);
      dxy = gxx(i,1);
      dyy = gxx(i,2);
      for (int k=0; k<nqp; k++){
	k3           = k2+k;
	pxx[k1x +k] += p[k3]*dx;
	pxx[k1y +k] += p[k3]*dy;
	pxx[k1xx+k] += p[k3]*dxx;
	pxx[k1xy+k] += p[k3]*dxy;
	pxx[k1yy+k] += p[k3]*dyy;
      }}}
}
