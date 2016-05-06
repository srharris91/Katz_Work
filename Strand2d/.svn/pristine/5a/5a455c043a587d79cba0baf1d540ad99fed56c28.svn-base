#include "StrandBlockSolver.h"


void StrandBlockSolver::coarseMetrics()
{
  for (int n=0; n<nFaces; n++)
    for (int j=0; j<nPstr+1; j++)
      for (int m=0; m<ndim; m++) facu(m,j,n) = 0.;
  for (int n=0; n<nEdges; n++)
    for (int j=0; j<nPstr+1; j++)
      for (int m=0; m<ndim; m++) facs(m,j,n) = 0.;
  for (int n=0; n<nFaces+nBedges; n++)
    for (int j=0; j<nPstr+2; j++) v(j,n) = 0.;

  int i,jp,m,im,k,l;
  double ir;
  for (int n=0; n<nFacesP-nGfacesP; n++){
    i                = f2cc(n);
    facu(0,0    ,i) += (*facuP)(0,0     ,n);
    facu(1,0    ,i) += (*facuP)(1,0     ,n);
    facu(0,nPstr,i) += (*facuP)(0,nPstrP,n);
    facu(1,nPstr,i) += (*facuP)(1,nPstrP,n);
    k                = 1;
    for (int j=1; j<nPstrP; j++){
    jp           = j+1;
  if (f2cs(j) != f2cs(jp)){
    facu(0,k,i) += (*facuP)(0,j,n);
    facu(1,k,i) += (*facuP)(1,j,n);
    k++;
  }}}

  for (int n=0; n<nEdgesP; n++){
    i  = f2ce(n);
    ir = 1.;
    if (i < 0) ir = -1.;
    i  = fabs(i);
  if (i != 0){
    im = i-1; // f2ce is 1-based to allow for meaningful signs
  for (int j=1; j<nPstrP+1; j++){
    k = f2cs(j);
    facs(0,k,im) += ir*(*facsP)(0,j,n);
    facs(1,k,im) += ir*(*facsP)(1,j,n);
  }}}

  k = nPstr+1;
  m = nPstrP+1;
  for (int n=0; n<nFacesP+nBedgesP; n++){
    i       = f2cc(n);
  if (i >= 0){
    v(0,i) += (*vP)(0,n);
    v(k,i) += (*vP)(m,n);
  for (int j=1; j<nPstrP+1; j++){
    l       = f2cs(j);
    v(l,i) += (*vP)(j,n);
  }}}
}
