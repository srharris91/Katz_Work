#include "Strand2dFCBlockSolver.h"
#include "solutionPoints1D.h"
#include "lagrangePoly1D.h"


void Strand2dFCBlockSolver::mapping()
{
  // set spacing in the mapping
  deltaN = 1./(double)(nStrandNode-1);
  deltaS = 2./(double)meshOrder;


  // allocate mapping terms
  xn.allocate(nSurfNode,nStrandNode);
  yn.allocate(nSurfNode,nStrandNode);
  xs.allocate(nSurfElem,meshOrder+1,nStrandNode);
  ys.allocate(nSurfElem,meshOrder+1,nStrandNode);
  jac.allocate(nSurfElem,meshOrder+1,nStrandNode);
  v.allocate(nSurfNode,nStrandNode);


  // determine mapping terms in the strand direction
  xn.set(0.);
  yn.set(0.);

  // strand mapping terms
  int n=0,j1,nj; //note: choose n=0 and then copy to all strands
  for (int j=0; j<nStrandNode; j++){
    j1 = icn2[j][0];
    nj = icn2[j][1];
    for (int m=0; m<nj; m++){
      xn(n,j) += icn1[j][m]*strandX(j1);
      j1++;
    }}
  for (int n=0; n<nSurfNode; n++)
    for (int j=0; j<nStrandNode; j++){
      xn(n,j) = xn(0,j);
      yn(n,j) = xn(0,j);
    }
  double nx,ny,dnr=1./deltaN;
  for (int n=0; n<nSurfNode; n++){
    nx = pointingVec(n,0);
    ny = pointingVec(n,1);
    for (int j=0; j<nStrandNode; j++){
      xn(n,j) *= (nx*dnr);
      yn(n,j) *= (ny*dnr);
    }}


  // determine mapping terms in the surface direction
  xs.set(0.);
  ys.set(0.);
  int spacing=0; // assume equally spaced points in surface elements for now
  Array1D<double> ss(meshOrder+1);
  solutionPoints1D(meshOrder,
		   spacing,
		   &ss(0));

  bool test=true;
  Array2D<double> lc(meshOrder+1,meshOrder+1);
  lagrangePoly1D(test, // coefficients to form Lagrange polynomials
		 meshOrder,
		 &ss(0),
		 &lc(0,0));

  // ls(i,j) = (dl_j/ds)_i (a row is all Lagrange polynomials (derivatives)
  // evaluated at a single mesh point i)
  ls.allocate(meshOrder+1,meshOrder+1);
  ls.set(0.);
  int km;
  for (int i=0; i<meshOrder+1; i++) // ith mesh point
    for (int j=0; j<meshOrder+1; j++) // jth Lagrange polynomial
      for (int k=0; k<meshOrder+1; k++){
	km       = max(0,k-1);
	ls(i,j) +=((double)k)*pow(ss(i),km)*lc(j,k);
      }

  /*
  // exact mapping for a circle (used for testing)
  int nn;
  double pi=4*atan(1.),t,dt=2.*pi/(double)nSurfElem;
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<meshOrder+1; i++){
      nn = surfElem(n,i);
      t  = atan2(surfX(nn,1),surfX(nn,0));
      for (int j=0; j<nStrandNode; j++){
	xs(n,i,j) = (.5+strandX(j))*sin(t)*.5*dt;
	ys(n,i,j) =-(.5+strandX(j))*cos(t)*.5*dt;
      }}
  */

  xs.set(0.);
  ys.set(0.);
  int ni,nm;
  double x0,y0;
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<meshOrder+1; i++){ // ith point in the element
      ni = surfElem(n,i);
      for (int m=0; m<meshOrder+1; m++){ // mth Lagrange poly. in mapping
	nm = surfElem(n,m);
	x0 = surfX(nm,0);
	y0 = surfX(nm,1);
	nx = pointingVec(nm,0);
	ny = pointingVec(nm,1);
	for (int j=0; j<nStrandNode; j++){
	  xs(n,i,j) += ls(i,m)*(x0+nx*strandX(j));
	  ys(n,i,j) += ls(i,m)*(y0+ny*strandX(j));
	}}
      for (int j=0; j<nStrandNode; j++)
	jac(n,i,j) = xs(n,i,j)*yn(ni,j)-ys(n,i,j)*xn(ni,j);
    }


  // compute nodal "volumes," which are just the sum of the Jacobians
  // at collocated nodes
  int i1,i2,n1,n2;
  v.set(0.);
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<nElemEdge; i++){
      i1 = elemEdge(i,0);
      i2 = elemEdge(i,1);
      n1 = surfElem(n,i1);
      n2 = surfElem(n,i2);
      for (int j=0; j<nStrandNode; j++){
	v(n1,j) += .5*deltaS*jac(n,i1,j);
	v(n2,j) += .5*deltaS*jac(n,i2,j);
      }}


  // compute agglomerated metric terms for the flux in the strand direction
  xsA.allocate(nSurfNode,nStrandNode);
  ysA.allocate(nSurfNode,nStrandNode);
  xsA.set(0.);
  ysA.set(0.);
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<nElemEdge; i++){
      i1 = elemEdge(i,0);
      i2 = elemEdge(i,1);
      n1 = surfElem(n,i1);
      n2 = surfElem(n,i2);
      for (int j=0; j<nStrandNode; j++){
	xsA(n1,j) += .5*deltaS*xs(n,i1,j);
	ysA(n1,j) += .5*deltaS*ys(n,i1,j);
	xsA(n2,j) += .5*deltaS*xs(n,i2,j);
	ysA(n2,j) += .5*deltaS*ys(n,i2,j);
      }}


  // s-derivatives of xn and yn (needed for flux gradients
  xns.allocate(nSurfElem,meshOrder+1,nStrandNode);
  yns.allocate(nSurfElem,meshOrder+1,nStrandNode);
  xns.set(0.);
  yns.set(0.);
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<meshOrder+1; i++) // ith point in the element
      for (int m=0; m<meshOrder+1; m++){ // mth Lagrange poly. in mapping
	nm = surfElem(n,m);
	for (int j=0; j<nStrandNode; j++){
	  xns(n,i,j) += ls(i,m)*xn(nm,j);
	  yns(n,i,j) += ls(i,m)*yn(nm,j);
	}}


  // clean up
  ss.deallocate();
  lc.deallocate();
}
