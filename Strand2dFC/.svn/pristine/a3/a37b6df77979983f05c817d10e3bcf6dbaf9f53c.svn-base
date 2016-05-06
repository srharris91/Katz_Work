#include "Strand2dFCBlockSolver.h"
#include "lagrangePoly1D.h"
#include "solutionPoints1D.h"


void Strand2dFCBlockSolver::initRestrict(const int& meshOrder0F0,
					 int* surfElem0F0,
					 int* bndNodeF0,
					 double* qF0,
					 double* rF0)
{
  // copy needed data from the fine level
  meshOrder0F = meshOrder0F0;
  surfElem0F  = surfElem0F0;
  bndNodeF    = bndNodeF0;
  qF          = qF0;
  rF          = rF0;


  // form multigrid solution restriction stencils
  psp2MGS.allocate(nSurfNode+1);
  psp2MGS.set(0);
  for (int n=0; n<nSurfElem0; n++){
    psp2MGS(surfElem0(n,0)+1) = 1;
    psp2MGS(surfElem0(n,1)+1) = 1;
    for (int i=2; i<meshOrder0+1; i++)
      psp2MGS(surfElem0(n,i)+1) = meshOrder0F+1;
  }

  for(int n=1; n<nSurfNode+1; n++) psp2MGS(n) += psp2MGS(n-1);
  npsp1MGS = psp2MGS(nSurfNode);
  psp1MGS.allocate(npsp1MGS);
  wsp1MGS.allocate(npsp1MGS);

  int m;
  Array1D<int> flag(nSurfNode);
  flag.set(-1);
  for (int n=0; n<nSurfElem0; n++){
    m = surfElem0(n,0);
    if (flag(m) == -1){
      psp1MGS(psp2MGS(m)++) = surfElem0F[n*(meshOrder0F+1)  ];
      flag(m)               = 0;
    }
    m = surfElem0(n,1);
    if (flag(m) == -1){
      psp1MGS(psp2MGS(m)++) = surfElem0F[n*(meshOrder0F+1)+1];
      flag(m)               = 0;
    }
    for (int i=2; i<meshOrder0+1; i++){
      m = surfElem0(n,i);
      for (int ii=0; ii<meshOrder0F+1; ii++)
	psp1MGS(psp2MGS(m)++) = surfElem0F[n*(meshOrder0F+1)+ii];
    }}

  for(int n=nSurfNode; n>0; n--) psp2MGS(n) = psp2MGS(n-1);
  psp2MGS(0) = 0;


  // form solution restriction weights, which are interpolation coefficients
  int spacing=0; // assume equally spaced points in surface elements for now
  Array1D<double> ss(meshOrder0+1);
  solutionPoints1D(meshOrder0,
		   spacing,
		   &ss(0));
  Array1D<double> ssF(meshOrder0F+1);
  solutionPoints1D(meshOrder0F,
		   spacing,
		   &ssF(0));
  bool test=false;
  Array2D<double> lcF(meshOrder0F+1,meshOrder0F+1);
  lagrangePoly1D(test, // coefficients to form Lagrange polynomials
		 meshOrder0F,
		 &ssF(0),
		 &lcF(0,0));

  // l(i,j) = (l_j)_i (a row is all fine level Lagrange polynomials
  // evaluated at a single coarse level mesh point i)
  Array2D<double> l(meshOrder0+1,meshOrder0F+1);
  l.set(0.);
  for (int i=0; i<meshOrder0+1; i++) // ith mesh point
    for (int j=0; j<meshOrder0F+1; j++) // jth Lagrange polynomial
      for (int k=0; k<meshOrder0F+1; k++)
	l(i,j) += pow(ss(i),k)*lcF(j,k);

  int ni,nSurfNodeF=nSurfElem0*meshOrder0F+nBndNode/2;
  Array1D<double> a(nSurfNodeF);
  a.set(0.);
  for (int n=0; n<nSurfElem0; n++){
    m = surfElem0(n,0);
    for (int ii=0; ii<meshOrder0F+1; ii++){
      ni    = surfElem0F[n*(meshOrder0F+1)+ii];
      a(ni) = l(0,ii);
    }
    for (int ii=psp2MGS(m); ii<psp2MGS(m+1); ii++)
      wsp1MGS(ii) = a(psp1MGS(ii));
    for (int ii=0; ii<meshOrder0F+1; ii++){
      ni    = surfElem0F[n*(meshOrder0F+1)+ii];
      a(ni) = 0.;
    }

    m = surfElem0(n,1);
    for (int ii=0; ii<meshOrder0F+1; ii++){
      ni    = surfElem0F[n*(meshOrder0F+1)+ii];
      a(ni) = l(1,ii);
    }
    for (int ii=psp2MGS(m); ii<psp2MGS(m+1); ii++)
      wsp1MGS(ii) = a(psp1MGS(ii));
    for (int ii=0; ii<meshOrder0F+1; ii++){
      ni    = surfElem0F[n*(meshOrder0F+1)+ii];
      a(ni) = 0.;
    }

    for (int i=2; i<meshOrder0+1; i++){
      m = surfElem0(n,i);
      for (int ii=0; ii<meshOrder0F+1; ii++){
	ni    = surfElem0F[n*(meshOrder0F+1)+ii];
	a(ni) = l(i,ii);
      }
      for (int ii=psp2MGS(m); ii<psp2MGS(m+1); ii++)
	wsp1MGS(ii) = a(psp1MGS(ii));
      for (int ii=0; ii<meshOrder0F+1; ii++){
	ni    = surfElem0F[n*(meshOrder0F+1)+ii];
	a(ni) = 0.;
      }}}

  /*
  for (int n=0; n<nSurfNode; n++){
    cout << "\nSurfNode: " << n << " " << psp2MGS(n+1)-psp2MGS(n) << endl;
    for (int i=psp2MGS(n); i<psp2MGS(n+1); i++)
      cout << psp1MGS(i) << " " << wsp1MGS(i) << endl;
  }
  if (level == 2) exit(0);
  */


  // form multigrid residual restriction stencils
  psp2MGR.allocate(nSurfNodeF+1);
  psp2MGR.set(0);
  for (int n=0; n<nSurfElem0; n++){
    m            = surfElem0F[n*(meshOrder0F+1)  ];
    psp2MGR(m+1) = 1;
    m            = surfElem0F[n*(meshOrder0F+1)+1];
    psp2MGR(m+1) = 1;
    for (int i=2; i<meshOrder0F+1; i++){
      m            = surfElem0F[n*(meshOrder0F+1)+i];
      psp2MGR(m+1) = meshOrder0+1;
    }}

  for(int n=1; n<nSurfNodeF+1; n++) psp2MGR(n) += psp2MGR(n-1);
  npsp1MGR = psp2MGR(nSurfNodeF);
  psp1MGR.allocate(npsp1MGR);
  wsp1MGR.allocate(npsp1MGR);

  flag.deallocate();
  flag.allocate(nSurfNodeF);
  flag.set(-1);
  for (int n=0; n<nSurfElem0; n++){
    m = surfElem0F[n*(meshOrder0F+1)  ];
    if (flag(m) == -1){
      psp1MGR(psp2MGR(m)++) = surfElem0(n,0);
      flag(m)               = 0;
    }
    m = surfElem0F[n*(meshOrder0F+1)+1];
    if (flag(m) == -1){
      psp1MGR(psp2MGR(m)++) = surfElem0(n,1);
      flag(m)               = 0;
    }
    for (int i=2; i<meshOrder0F+1; i++){
      m = surfElem0F[n*(meshOrder0F+1)+i];
      for (int ii=0; ii<meshOrder0+1; ii++)
	psp1MGR(psp2MGR(m)++) = surfElem0(n,ii);
    }}

  for(int n=nSurfNodeF; n>0; n--) psp2MGR(n) = psp2MGR(n-1);
  psp2MGR(0) = 0;


  // form residual restriction weights, which are interpolation coefficients
  Array2D<double> lc(meshOrder0+1,meshOrder0+1);
  lagrangePoly1D(test, // coefficients to form Lagrange polynomials
		 meshOrder0,
		 &ss(0),
		 &lc(0,0));

  // l(i,j) = (l_j)_i (a row is all coarse level Lagrange polynomials
  // evaluated at a single fine level mesh point i)
  l.deallocate();
  l.allocate(meshOrder0F+1,meshOrder0+1);
  l.set(0.);
  for (int i=0; i<meshOrder0F+1; i++) // ith mesh point
    for (int j=0; j<meshOrder0+1; j++) // jth Lagrange polynomial
      for (int k=0; k<meshOrder0+1; k++)
	l(i,j) += pow(ssF(i),k)*lc(j,k);

  a.deallocate();
  a.allocate(nSurfNode);
  a.set(0.);
  for (int n=0; n<nSurfElem0; n++){
    m = surfElem0F[n*(meshOrder0F+1)  ];
    for (int ii=0; ii<meshOrder0+1; ii++){
      ni    = surfElem0(n,ii);
      a(ni) = l(0,ii);
    }
    for (int ii=psp2MGR(m); ii<psp2MGR(m+1); ii++)
      wsp1MGR(ii) = a(psp1MGR(ii));
    for (int ii=0; ii<meshOrder0+1; ii++){
      ni    = surfElem0(n,ii);
      a(ni) = 0.;
    }

    m = surfElem0F[n*(meshOrder0F+1)+1];
    for (int ii=0; ii<meshOrder0+1; ii++){
      ni    = surfElem0(n,ii);
      a(ni) = l(1,ii);
    }
    for (int ii=psp2MGR(m); ii<psp2MGR(m+1); ii++)
      wsp1MGR(ii) = a(psp1MGR(ii));
    for (int ii=0; ii<meshOrder0+1; ii++){
      ni    = surfElem0(n,ii);
      a(ni) = 0.;
    }

    for (int i=2; i<meshOrder0F+1; i++){
      m = surfElem0F[n*(meshOrder0F+1)+i];
      for (int ii=0; ii<meshOrder0+1; ii++){
	ni    = surfElem0(n,ii);
	a(ni) = l(i,ii);
      }
      for (int ii=psp2MGR(m); ii<psp2MGR(m+1); ii++)
	wsp1MGR(ii) = a(psp1MGR(ii));
      for (int ii=0; ii<meshOrder0+1; ii++){
	ni    = surfElem0(n,ii);
	a(ni) = 0.;
      }}}


  /*
  // set contributions to boundary nodes to 0.
  for (int n=0; n<nSurfNodeF; n++)
    for (int i=psp2MGR(n); i<psp2MGR(n+1); i++){
      m = psp1MGR(i);
      for (int ii=0; ii<nBndNode; ii++)
	if (m == bndNode(ii)) wsp1MGR(i) = 0.;
    }
  */

  /*
  for (int n=0; n<nSurfNodeF; n++){
    cout << "\nSurfNode: " << n << " " << psp2MGR(n+1)-psp2MGR(n) << endl;
    for (int i=psp2MGR(n); i<psp2MGR(n+1); i++)
      cout << psp1MGR(i) << " " << wsp1MGR(i) << endl;
  }
  if (level == 1) exit(0);
  */


  // clean up
  flag.deallocate();
  ss.deallocate();
  ssF.deallocate();
  lcF.deallocate();
  lc.deallocate();
  l.deallocate();
  a.deallocate();
}
