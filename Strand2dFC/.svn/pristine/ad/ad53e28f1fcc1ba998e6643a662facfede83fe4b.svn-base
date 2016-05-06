#include "lagrangePoly.h"
#include "Array2D.h"
#include "matinv.h"
#include <math.h>
#include <iostream>

// determine coefficients for each Lagrange polynomial
// lc(i,j) is the jth coefficient of the ith Lagrange polynomial
// (a row is a complete Lagrange polynomial)


void lagrangePoly(const bool& test,
		  const int& order,
		  const double* rs,
		  double* lc)
{
  int j,m,nsp=(order+2)*(order+1)/2;
  Array2D<double> A(nsp,nsp);
  for (int n=0; n<nsp; n++){ // for each Lagrange polynomial
    for (int i=0; i<nsp; i++){ // evaluate at each solution point
      m = 3*i;
      j = 0;
      for (int k=0; k<=order; k++)
	for (int l=0; l<=order-k; l++)
	  A(i,j++) = pow(rs[m],k)*pow(rs[m+1],l);
    }
    matinv(nsp,&A(0,0)); // compute inv(A) to get Lagrange coefficients
    for (int i=0; i<nsp; i++) lc[n*nsp+i] = A(i,n);
  }


  // test Lagrange polynomial delta property
  if (test){
    std::cout << "\nTest for Lagrange polynomial delta property: " << std::endl;
    for (int n=0; n<nsp; n++){
      A.set(0.);
      for (int i=0; i<nsp; i++){ // evaluate at each solution point
	m = 3*i;
	j = 0;
	for (int k=0; k<=order; k++)
	  for (int l=0; l<=order-k; l++)
	    A(0,i) += pow(rs[m],k)*pow(rs[m+1],l)*lc[n*nsp+j++];
      }
      for (int i=0; i<nsp; i++){
	if (fabs(A(0,i)) < 5.e-13) A(0,i) = 0.;
	std::cout << A(0,i) << " ";
      }
      std::cout << std::endl;
    }}

  A.deallocate();
}
