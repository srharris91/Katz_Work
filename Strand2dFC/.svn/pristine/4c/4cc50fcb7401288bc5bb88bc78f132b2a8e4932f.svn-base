#include "Strand2dFCBlockSolver.h"
#include "solutionPoints1D.h"
#include "lagrangePoly1D.h"


void Strand2dFCBlockSolver::sourceSetup()
{
  // initialize coefficient array
  wsp1S.allocate(npsp1S);
  wsp1S.set(0.);
  int i1,i2,n1,n2;
  double ds5=.5*deltaS,ds8=.125*deltaS*deltaS,ds6=deltaS/6.;


  // use elements to form source coefficients
  Array2D<double> a(nSurfElem,meshOrder+1);
  a.set(0.);
  if (surfOrder == 1 || surfOrder == 2) //mass lumped
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<nElemEdge; i++){
	i1    = elemEdge(i,0);
	i2    = elemEdge(i,1);
	n1    = surfElem(n,i1);
	n2    = surfElem(n,i2);

	// contributions to n1
	a(n,i1) = .5*deltaS;
	for (int m=psp2S(n1); m<psp2S(n1+1); m++) //add to coefficient array
	  wsp1S(m) += a(psp1S(m,0),psp1S(m,1));
	a(n,i1) = 0.; //reset to 0

	// contributions to n2
	a(n,i2) = .5*deltaS;
	for (int m=psp2S(n2); m<psp2S(n2+1); m++) //add to coefficient array
	  wsp1S(m) += a(psp1S(m,0),psp1S(m,1));
	a(n,i2) = 0.; //reset to 0
      }

  /*
  else if (surfOrder == 2) //Galerkin
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<nElemEdge; i++){
	i1 = elemEdge(i,0);
	i2 = elemEdge(i,1);
	n1 = surfElem(n,i1);
	n2 = surfElem(n,i2);

	// contributions to n1
	a(n,i1) = ds6*2.;
	a(n,i2) = ds6;
	for (int m=psp2S(n1); m<psp2S(n1+1); m++) //add to coefficient array
	  wsp1S(m) += a(psp1S(m,0),psp1S(m,1));
	a(n,i1) = 0.; //reset to 0
	a(n,i2) = 0.;

	// contributions to n2
	a(n,i1) = ds6;
	a(n,i2) = ds6*2.;
	for (int m=psp2S(n2); m<psp2S(n2+1); m++) //add to coefficient array
	  wsp1S(m) += a(psp1S(m,0),psp1S(m,1));
	a(n,i1) = 0.; //reset to 0
	a(n,i2) = 0.;
      }
  */

  else if (surfOrder == 3) //corrected source by edge
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<nElemEdge; i++){
	i1 = elemEdge(i,0);
	i2 = elemEdge(i,1);
	n1 = surfElem(n,i1);
	n2 = surfElem(n,i2);

	// contributions to n1
	a(n,i1) = ds6*2.;
	a(n,i2) = ds6;
	for (int m=0; m<meshOrder+1; m++){
	  a(n,m) -= ds6*ds5*(ls(i1,m)+ls(i2,m));
	  for (int mm=0; mm<meshOrder+1; mm++)
	    a(n,mm) -= ds6*ds8*(ls(i1,m)*ls(m,mm)+ls(i2,m)*ls(m,mm));
	}
	for (int m=psp2S(n1); m<psp2S(n1+1); m++) //add to coefficient array
	  wsp1S(m) += a(psp1S(m,0),psp1S(m,1));
	for (int m=0; m<meshOrder+1; m++) a(n,m) = 0.; //reset to 0

	// contributions to n2
	a(n,i1) = ds6;
	a(n,i2) = ds6*2.;
	for (int m=0; m<meshOrder+1; m++){
	  a(n,m) += ds6*ds5*(ls(i1,m)+ls(i2,m));
	  for (int mm=0; mm<meshOrder+1; mm++)
	    a(n,mm) -= ds6*ds8*(ls(i1,m)*ls(m,mm)+ls(i2,m)*ls(m,mm));
	}
	for (int m=psp2S(n2); m<psp2S(n2+1); m++) //add to coefficient array
	  wsp1S(m) += a(psp1S(m,0),psp1S(m,1));
	for (int m=0; m<meshOrder+1; m++) a(n,m) = 0.; //reset to 0
      }


  // clean up
  a.deallocate();

  /*
  double sum;
  for (int n=0; n<nSurfNode; n++){
    cout << "\nSurfNode: " << n << " " << psp2S(n+1)-psp2S(n) << endl;
    sum = 0.;
    for (int i=psp2S(n); i<psp2S(n+1); i++){
      sum += wsp1S(i);
      cout << psp1S(i,0)  << " " << psp1S(i,1) << " "
	   << wsp1S(i) << endl;
    }
    cout << sum << endl;
  }
  exit(0);
  */
}
