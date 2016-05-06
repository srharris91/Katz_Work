#include "StrandBlockSolver.h"


void StrandBlockSolver::perturbNodes()
{
  // find min distance of surrounding nodes
  int jj=nPstr+1,n1,n2,jm;
  double dx,dy,ds;
  Array2D<double> dsmin(jj,nNodes);

  for (int n=0; n<nNodes; n++)
    for (int j=0; j<nPstr+1; j++) dsmin(j,n) = 1.e14;

  for (int n=0; n<nFaces; n++){
    n1 = face(0,n);
    n2 = face(1,n);
  for (int j=0; j<nPstr+1; j++){
    dx = x(0,j,n1)-x(0,j,n2);
    dy = x(1,j,n1)-x(1,j,n2);
    ds = dx*dx+dy*dy;
    if (ds < dsmin(j,n1)) dsmin(j,n1) = ds;
    if (ds < dsmin(j,n2)) dsmin(j,n2) = ds;
  }}

  for (int n=0; n<nNodes; n++){
  for (int j=1; j<nPstr+1; j++){
    jm = j-1;
    dx = x(0,j,n)-x(0,jm,n);
    dy = x(1,j,n)-x(1,jm,n);
    ds = dx*dx+dy*dy;
    if (ds < dsmin(j, n)) dsmin(j, n) = ds;
    if (ds < dsmin(jm,n)) dsmin(jm,n) = ds;
  }}

  for (int n=0; n<nNodes; n++)
    for (int j=0; j<nPstr+1; j++) dsmin(j,n) = sqrt(dsmin(j,n));


  // mark boundary nodes
  int flag[nNodes];
  for (int n=0; n<nNodes; n++) flag[n] = 0;
  for (int n=nEdges-nBedges; n<nEdges; n++){
    n1       = edgn(n);
    flag[n1] = 1;
  }


  // perturb nodes
  double pi = 4.*atan(1.),fact=.1;
  for (int n=0; n<nNodes-nGnodes; n++){
    if (flag[n] == 0){
      for (int j=1; j<nPstr; j++){
	dx       = double(rand())/double(RAND_MAX);
	dy       = double(rand())/double(RAND_MAX);
	dx       = fact*dsmin(j,n)*dx;
	dy       = 2.*pi*dy;
	x(0,j,n) = x(0,j,n)+dx*cos(dy);
	x(1,j,n) = x(1,j,n)+dx*sin(dy);
      }}}

  dsmin.deallocate();
}
