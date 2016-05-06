#include "Tri2dFCBlockSolver.h"
#include "lagrangePoly.h"


void Tri2dFCBlockSolver::initProlong(int levelC0,
				     int* elemC0,
				     double* q0C0,
				     double* qC0)
{
  // set the pointers to the coarse data
  levelC = levelC0;
  elemC  = elemC0;
  q0C    = q0C0;
  qC     = qC0;


  // for each fine node, find the coarse level element which contains it,
  // as well as its (r,s) index in the standard element
  nce.allocate(nNode,2);
  nce.set(-1);
  for (int n=0; n<nElem; n++)
    for (int j=0; j<nne; j++)
      if (nce(elem(n,j),0) == -1){
	nce(elem(n,j),0) = n;
	nce(elem(n,j),1) = j;
      }


  // set up solution restriction coefficients via interpolation
  // lqCF(i,j) = l_j(r_i) (a row is all coarse level Lagrange polynomials
  // evaluated at a single fine level solution point i)
  int orderE=3-level; //order of local elements on this level
  int orderEC=3-levelC; //order of local elements on coarse level
  bool test=false;
  nneC =(orderEC+2)*(orderEC+1)/2; //number of nodes per element on coarse level
  Array2D<double> rs(nne,3),rsC(nneC,3),lcC(nneC,nneC);
  solutionPoints(orderE,
		 spacing,
		 &rs(0,0));
  solutionPoints(orderEC,
		 spacing,
		 &rsC(0,0));
  lagrangePoly(test,
	       orderEC,
	       &rsC(0,0),
	       &lcC(0,0));
  lqCF.allocate(nne,nneC);
  lqCF.set(0.);
  int j;
  double ri,si;
  for (int n=0; n<nneC; n++) // nth coarse level Lagrange polynomial
    for (int i=0; i<nne; i++){ // ith coarse level solution point
      j = 0;
      for (int k=0; k<=orderEC; k++)
	for (int l=0; l<=orderEC-k; l++){
	  ri = rs(i,0);
	  si = rs(i,1);
	  lqCF(i,n) += pow(ri,k)*pow(si,l)*lcC(n,j++);
	}
    }

  rs.deallocate();
  rsC.deallocate();
  lcC.deallocate();

  /*
  for (int i=0; i<nne; i++){
    for (int j=0; j<nneC; j++) cout << lqCF(i,j) << " ";
    cout << endl;
  }
  */
}
