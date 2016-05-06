#include "Strand2dFCBlockSolver.h"
#include "solutionPoints1D.h"
#include "lagrangePoly1D.h"


void Strand2dFCBlockSolver::gradSetupFull(Array3D<double>& gx)
{
  // initialize coefficient array
  gx.set(0.);


  // degree at each surface node
  Array1D<int> sum(nSurfNode);
  sum.set(0);
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<meshOrder+1; i++) sum(surfElem(n,i))++;


  // use elements to form gradient coefficients
  int ni,nm;
  double xnj,ynj;
  Array3D<double> a(nSurfNode,nStrandNode,2);
  a.set(0.);
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<meshOrder+1; i++){ //ith point in the element
      ni = surfElem(n,i);
      for (int m=0; m<meshOrder+1; m++){ //mth Lagrange poly. in mapping
	nm = surfElem(n,m);
	for (int j=0; j<nStrandNode; j++){ //local gradient coefficients
	  xnj       = xn(ni,j)/(jac(n,i,j)*(double)sum(ni));
	  ynj       = yn(ni,j)/(jac(n,i,j)*(double)sum(ni));
	  a(nm,j,0) = ls(i,m)*ynj;
	  a(nm,j,1) =-ls(i,m)*xnj;
	}}
      for (int m=psp2(ni); m<psp2(ni+1); m++){ //add to coefficient array
	nm = psp1(m);
	for (int j=0; j<nStrandNode; j++){
	  gx(m,j,0) += a(nm,j,0);
	  gx(m,j,1) += a(nm,j,1);
	}}
      for (int m=0; m<meshOrder+1; m++){ //reset helper array to zero
	nm = surfElem(n,m);
	for (int j=0; j<nStrandNode; j++){
	  a(nm,j,0) = 0.;
	  a(nm,j,1) = 0.;
	}}}


  // clean up
  sum.deallocate();
  a.deallocate();
}
