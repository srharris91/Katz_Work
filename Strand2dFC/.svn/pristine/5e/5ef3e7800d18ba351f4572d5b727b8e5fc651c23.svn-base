#include "Strand2dFCBlockSolver.h"
#include "lagrangePoly1D.h"
#include "solutionPoints1D.h"


void Strand2dFCBlockSolver::initProlong(const int& meshOrder0C0,
					int* surfElem0C0,
					int* bndNodeC0,
					double* q0C0,
					double* qC0)
{
  // copy needed data from the coarse level
  meshOrder0C = meshOrder0C0;
  surfElem0C  = surfElem0C0;
  bndNodeC    = bndNodeC0;
  q0C         = q0C0;
  qC          = qC0;


  // form multigrid correction prolongation stencils
  psp2MGC.allocate(nSurfNode+1);
  psp2MGC.set(0);
  for (int n=0; n<nSurfElem0; n++){
    psp2MGC(surfElem0(n,0)+1) = 1;
    psp2MGC(surfElem0(n,1)+1) = 1;
    for (int i=2; i<meshOrder0+1; i++)
      psp2MGC(surfElem0(n,i)+1) = meshOrder0C+1;
  }

  for(int n=1; n<nSurfNode+1; n++) psp2MGC(n) += psp2MGC(n-1);
  npsp1MGC = psp2MGC(nSurfNode);
  psp1MGC.allocate(npsp1MGC);
  wsp1MGC.allocate(npsp1MGC);

  int m;
  Array1D<int> flag(nSurfNode);
  flag.set(-1);
  for (int n=0; n<nSurfElem0; n++){
    m = surfElem0(n,0);
    if (flag(m) == -1){
      psp1MGC(psp2MGC(m)++) = surfElem0C[n*(meshOrder0C+1)  ];
      flag(m)               = 0;
    }
    m = surfElem0(n,1);
    if (flag(m) == -1){
      psp1MGC(psp2MGC(m)++) = surfElem0C[n*(meshOrder0C+1)+1];
      flag(m)               = 0;
    }
    for (int i=2; i<meshOrder0+1; i++){
      m = surfElem0(n,i);
      for (int ii=0; ii<meshOrder0C+1; ii++)
	psp1MGC(psp2MGC(m)++) = surfElem0C[n*(meshOrder0C+1)+ii];
    }}

  for(int n=nSurfNode; n>0; n--) psp2MGC(n) = psp2MGC(n-1);
  psp2MGC(0) = 0;


  // form residual restriction weights, which are interpolation coefficients
  int spacing=0; // assume equally spaced points in surface elements for now
  Array1D<double> ss(meshOrder0+1);
  solutionPoints1D(meshOrder0,
		   spacing,
		   &ss(0));
  Array1D<double> ssC(meshOrder0C+1);
  solutionPoints1D(meshOrder0C,
		   spacing,
		   &ssC(0));
  bool test=false;
  Array2D<double> lcC(meshOrder0C+1,meshOrder0C+1);
  lagrangePoly1D(test, // coefficients to form Lagrange polynomials
		 meshOrder0C,
		 &ssC(0),
		 &lcC(0,0));

  // l(i,j) = (l_j)_i (a row is all coarse level Lagrange polynomials
  // evaluated at a single fine level mesh point i)
  Array2D<double> l(meshOrder0+1,meshOrder0C+1);
  l.set(0.);
  for (int i=0; i<meshOrder0+1; i++) // ith mesh point
    for (int j=0; j<meshOrder0C+1; j++) // jth Lagrange polynomial
      for (int k=0; k<meshOrder0C+1; k++)
	l(i,j) += pow(ss(i),k)*lcC(j,k);

  int ni,nSurfNodeC=nSurfElem0*meshOrder0C+nBndNode/2;
  Array1D<double> a(nSurfNodeC);
  a.set(0.);
  for (int n=0; n<nSurfElem0; n++){
    m = surfElem0(n,0);
    for (int ii=0; ii<meshOrder0C+1; ii++){
      ni    = surfElem0C[n*(meshOrder0C+1)+ii];
      a(ni) = l(0,ii);
    }
    for (int ii=psp2MGC(m); ii<psp2MGC(m+1); ii++)
      wsp1MGC(ii) = a(psp1MGC(ii));
    for (int ii=0; ii<meshOrder0C+1; ii++){
      ni    = surfElem0C[n*(meshOrder0C+1)+ii];
      a(ni) = 0.;
    }

    m = surfElem0(n,1);
    for (int ii=0; ii<meshOrder0C+1; ii++){
      ni    = surfElem0C[n*(meshOrder0C+1)+ii];
      a(ni) = l(1,ii);
    }
    for (int ii=psp2MGC(m); ii<psp2MGC(m+1); ii++)
      wsp1MGC(ii) = a(psp1MGC(ii));
    for (int ii=0; ii<meshOrder0C+1; ii++){
      ni    = surfElem0C[n*(meshOrder0C+1)+ii];
      a(ni) = 0.;
    }

    for (int i=2; i<meshOrder0+1; i++){
      m = surfElem0(n,i);
      for (int ii=0; ii<meshOrder0C+1; ii++){
	ni    = surfElem0C[n*(meshOrder0C+1)+ii];
	a(ni) = l(i,ii);
      }
      for (int ii=psp2MGC(m); ii<psp2MGC(m+1); ii++)
	wsp1MGC(ii) = a(psp1MGC(ii));
      for (int ii=0; ii<meshOrder0C+1; ii++){
	ni    = surfElem0C[n*(meshOrder0C+1)+ii];
	a(ni) = 0.;
      }}}

  /*
  for (int n=0; n<nSurfNode; n++){
    cout << "\nSurfNode: " << n << " " << psp2MGC(n+1)-psp2MGC(n) << endl;
    for (int i=psp2MGC(n); i<psp2MGC(n+1); i++)
      cout << psp1MGC(i) << " " << wsp1MGC(i) << endl;
  }
  if (level == 1) exit(0);
  */


  // clean up
  flag.deallocate();
  ss.deallocate();
  ssC.deallocate();
  lcC.deallocate();
  l.deallocate();
  a.deallocate();
}
