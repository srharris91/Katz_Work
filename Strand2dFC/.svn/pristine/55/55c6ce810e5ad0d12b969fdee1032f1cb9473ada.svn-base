#include "Strand2dFCBlockSolver.h"
#include "solutionPoints1D.h"
#include "lagrangePoly1D.h"
#include "gauss1d.h"


void Strand2dFCBlockSolver::bNormal()
{
  // allocate and initialize surface normals (includes at strand tips)
  sn.allocate(nSurfNode,2,2);
  sn.set(0.);


  // compute surface and tip node normals
  int m;
  double xsi,ysi,ds;
  for (int n=0; n<nSurfElem; n++){
    for (int i=0; i<meshOrder+1; i++){ // strand roots
      m          = surfElem(n,i);
      xsi        = xs(n,i,0);
      ysi        = ys(n,i,0);
      ds         = 1./sqrt(xsi*xsi+ysi*ysi);
      xsi       *= ds;
      ysi       *= ds;
      sn(m,0,0) += ysi;
      sn(m,0,1) -= xsi;
    }
    for (int i=0; i<meshOrder+1; i++){ // strand tips
      m          = surfElem(n,i);
      xsi        = xs(n,i,nStrandNode-1);
      ysi        = ys(n,i,nStrandNode-1);
      ds         = 1./sqrt(xsi*xsi+ysi*ysi);
      xsi       *= ds;
      ysi       *= ds;
      sn(m,1,0) -= ysi;
      sn(m,1,1) += xsi;
    }}
  for (int n=0; n<nSurfNode; n++){
    ds         = 1./sqrt(sn(n,0,0)*sn(n,0,0)+sn(n,0,1)*sn(n,0,1));
    sn(n,0,0) *= ds;
    sn(n,0,1) *= ds;
    ds         = 1./sqrt(sn(n,1,0)*sn(n,1,0)+sn(n,1,1)*sn(n,1,1));
    sn(n,1,0) *= ds;
    sn(n,1,1) *= ds;
    /*
    cout << n
	 << " root: " << sn(n,0,0) << " " << sn(n,0,1) << " "
	 << " tip:  " << sn(n,1,0) << " " << sn(n,1,1) << endl;
    */
  }


  // form data structures for surface force integration
  // allocate data structures
  nQuadPoint = meshOrder;
  wQ.allocate(nQuadPoint);
  lQ.allocate(nQuadPoint,meshOrder+1);
  lsQ.allocate(nQuadPoint,meshOrder+1);
  xsQ.allocate(nSurfElem,nQuadPoint);
  ysQ.allocate(nSurfElem,nQuadPoint);
  xnQ.allocate(nSurfElem,nQuadPoint);
  ynQ.allocate(nSurfElem,nQuadPoint);
  jacQ.allocate(nSurfElem,nQuadPoint);

  Array2D<double> a(nQuadPoint,2);
  Array1D<double> sQ(nQuadPoint);
  gauss1d(meshOrder-1,
	  &a(0,0));
  for (int n=0; n<nQuadPoint; n++){
    sQ(n) = a(n,0);
    wQ(n) = a(n,1);
  }

  int spacing=0; // assume equally spaced points in surface elements for now
  Array1D<double> ss(meshOrder+1);
  solutionPoints1D(meshOrder,
		   spacing,
		   &ss(0));

  bool test=false;
  Array2D<double> lc(meshOrder+1,meshOrder+1);
  lagrangePoly1D(test, // coefficients to form Lagrange polynomials
		 meshOrder,
		 &ss(0),
		 &lc(0,0));

  // ls(i,j) = (dl_j/ds)_i (a row is all Lagrange polynomials (and derivatives)
  // evaluated at a single quadrature point i)
  lQ.set(0.);
  lsQ.set(0.);
  int km;
  for (int i=0; i<nQuadPoint; i++) // ith quadrature point
    for (int j=0; j<meshOrder+1; j++) // jth Lagrange polynomial
      for (int k=0; k<meshOrder+1; k++){
	km        = max(0,k-1);
	lQ(i,j)  +=            pow(sQ(i),k )*lc(j,k);
	lsQ(i,j) +=((double)k)*pow(sQ(i),km)*lc(j,k);
      }

  int nm;
  xnQ.set(0.);
  ynQ.set(0.);
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<nQuadPoint; i++)
      for (int m=0; m<meshOrder+1; m++){ // mth Lagrange poly. in mapping
	nm        = surfElem(n,m);
	xnQ(n,i) += lQ(i,m)*xn(nm,0);
	ynQ(n,i) += lQ(i,m)*yn(nm,0);
      }

  xsQ.set(0.);
  ysQ.set(0.);
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<nQuadPoint; i++){ // ith quadrature point
      for (int m=0; m<meshOrder+1; m++){ // mth Lagrange poly. in mapping
	nm        = surfElem(n,m);
	xsQ(n,i) += lsQ(i,m)*surfX(nm,0);
	ysQ(n,i) += lsQ(i,m)*surfX(nm,1);
      }
      jacQ(n,i) = xsQ(n,i)*ynQ(n,i)-ysQ(n,i)*xnQ(n,i);
    }


  // clean up
  a.deallocate();
  sQ.deallocate();
  ss.deallocate();
  lc.deallocate();
}
