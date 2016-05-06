#include "StrandBlockSolver.h"


void StrandBlockSolver::lhsTime(const int& step,
				const int& pseudoStep)
{
  // compute inviscid and viscous spectral radius
  specRadi();
  specRadv();


  // add physical time derivative to LHS
  if (step > 0){
    int jj=1,m;
    double tj[nq*nq],a=1.5/dtUnsteady,b;
    for (int n=0; n<nFaces-nGfaces; n++)
    for (int j=1; j<fClip(n)+1; j++){
      b = a*v(j,n);
      sys->lhsConsVarJacobian(jj,&q(0,j,n),&qa(0,j,n),&tj[0]);
      for (int k=0; k<nq; k++){
	m = k*nq;
      for (int l=0; l<nq; l++) dd(l,k,j,n) = tj[m+l]*b;
      }}}


  // compute the pseudo time step using sum of face areas around a cell
  // in place of volume, and add to LHS
  int j=nPstr+2,k=nFaces+nBedges,c1,c2,m;
  double Ax,Ay,A;
  Array2D<double> ll(j,k);
  for (int n=0; n<nFaces+nBedges; n++)
  for (int j=0; j<nPstr+2; j++) ll(j,n) = 0.;

  for (int n=0; n<nEdges; n++){
    c1       = edge(0,n);
    c2       = edge(1,n);
    m        = fClip(c1);
    if (fClip(c2) > m) m = fClip(c2);
  for (int j=1; j<m+1; j++){
    Ax       = facs(0,j,n);
    Ay       = facs(1,j,n);
    A        = Ax*Ax+Ay*Ay;
    ll(j,c1)+= A;
    ll(j,c2)+= A;
  }}
  int jp;
  for (int n=0; n<nFaces-nGfaces; n++){
  for (int j=0; j<fClip(n)+1; j++){
    jp       = j+1;
    Ax       = facu(0,j,n);
    Ay       = facu(1,j,n);
    A        = Ax*Ax+Ay*Ay;
    ll(j ,n)+= A;
    ll(jp,n)+= A;
  }}


  // CFL and VNN ramping
  double cflT,vnnT;
  if (pseudoStep > nRamp){
    cflT = cfl;
    vnnT = vnn;
  }
  else{
    cflT = cfl0+(cfl-cfl0)/((double)nRamp)*((double)pseudoStep);
    vnnT = vnn0+(vnn-vnn0)/((double)nRamp)*((double)pseudoStep);
  }

  for (int n=0; n<nFaces+nBedges; n++)
  for (int j=0; j<nPstr+2; j++) dt(j,n) = 0;
  if (inviscid == 1){
    for (int n=0; n<nFaces-nGfaces; n++)
    for (int j=0; j<fClip(n)+1; j++) dt(j,n) = radi(j,n)/cflT;
  }
  if (viscous == 1){
    double a;
    for (int n=0; n<nFaces-nGfaces; n++)
    for (int j=0; j<fClip(n)+1; j++){
      a = radv(j,n)/vnnT;
      if (a > dt(j,n)) dt(j,n) = a;
    }
  }
  for (int n=0; n<nFaces-nGfaces; n++)
    //for (int j=0; j<fClip(n)+1; j++) dt(j,n) = ll(j,n)/dt(j,n);
    for (int j=0; j<fClip(n)+1; j++) dt(j,n) = v(j,n)/dt(j,n);


  int jj=1;
  double pj[nq*nq],a;
  for (int n=0; n<nFaces-nGfaces; n++)
  for (int j=0; j<fClip(n)+1; j++){
    sys->lhsPreconJacobian(jj,&q(0,j,n),&qa(0,j,n),&pj[0]);
    a = v(j,n)/dt(j,n);
    for (int k=0; k<nq; k++){
      m = k*nq;
    for (int l=0; l<nq; l++) dd(l,k,j,n) += pj[m+l]*a;
    }}

  ll.deallocate();
}
