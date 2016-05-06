#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::gradient(const int& nqp,
				  const double* p,
				  double* px)
{
  for (int n=0; n<nNode*nqp*2; n++) px[n] = 0.;

  if (order == 2){ //Green-Gauss gradients
    int n1,n2,k1,k2,k1x,k1y,k2x,k2y,m=0;
    double Ax,Ay,a,third=1./3.;
    for (int n=0; n<nEdge; n++){ //interior edges
      n1  = edge(n,0);
      n2  = edge(n,1);
      Ax  = area(n,0);
      Ay  = area(n,1);
      k1  = n1*nqp;
      k2  = n2*nqp;
      k1x = n1*nqp*2;
      k1y = k1x+nqp;
      k2x = n2*nqp*2;
      k2y = k2x+nqp;
      for (int k=0; k<nqp; k++){
	a          = p[k1+k]+p[k2+k];
	px[k1x+k] += Ax*a;
	px[k1y+k] += Ay*a;
	px[k2x+k] -= Ax*a;
	px[k2y+k] -= Ay*a;
      }}
    for (int n=nEdge-nEdgeBd; n<nEdge; n++){ //boundary edges
      n1  = edge(n,0);
      n2  = edge(n,1);
      Ax  = areaBd(m  ,0);
      Ay  = areaBd(m++,1);
      k1  = n1*nqp;
      k2  = n2*nqp;
      k1x = n1*nqp*2;
      k1y = k1x+nqp;
      k2x = n2*nqp*2;
      k2y = k2x+nqp;
      for (int k=0; k<nqp; k++){
	a          =(5.*p[k1+k]+p[k2+k])*third;
	px[k1x+k] += Ax*a;
	px[k1y+k] += Ay*a;
	a          =(p[k1+k]+5.*p[k2+k])*third;
	px[k2x+k] += Ax*a;
	px[k2y+k] += Ay*a;
      }}
    for (int n=0; n<nNode; n++){ //normalize gradients
      a  = .5/v(n);
      k1 = n*nqp*2;
      for (int k=0; k<nqp*2; k++) px[k1+k] *= a;
    }}


  else if (order == 3){ //FEM gradients
    int k1x,k1y,k2,k3;
    double dx,dy;
    for (int n=0; n<nNode; n++){
      k1x = n*nqp*2;
      k1y = k1x+nqp;
      for(int i=psp2(n); i<psp2(n+1); i++){
	k2 = psp1(i)*nqp;
	dx = gx(i,0);
	dy = gx(i,1);
	for (int k=0; k<nqp; k++){
	  k3         = k2+k;
	  px[k1x+k] += p[k3]*dx;
	  px[k1y+k] += p[k3]*dy;
	}}}}
}
