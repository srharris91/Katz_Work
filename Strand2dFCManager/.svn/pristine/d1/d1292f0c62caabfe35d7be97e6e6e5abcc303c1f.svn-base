#include "lagrangePoly1D.h"
#include "Array2D.h"
#include "matinv.h"
#include <math.h>
#include <iostream>

// determine coefficients for each Lagrange polynomial
// lc(i,j) is the jth coefficient of the ith Lagrange polynomial
// (a row is a complete Lagrange polynomial)

void lagrangePoly1D(const bool& test,
		    const int& order,
		    const double* r,
		    double* lc)
{
  int j,m,nsp=order+1;
  Array2D<double> A(nsp,nsp);
  for (int n=0; n<nsp; n++){ // for each Lagrange polynomial
    for (int i=0; i<nsp; i++) // evaluate at each solution point
      for (int j=0; j<order+1; j++) A(i,j) = pow(r[i],j);
    matinv(nsp,&A(0,0)); // compute inv(A) to get Lagrange coefficients
    for (int i=0; i<nsp; i++) lc[n*nsp+i] = A(i,n);
  }


  // test Lagrange polynomial delta property
  if (test){
    std::cout << "\nTest for Lagrange polynomial delta property: " << std::endl;
    for (int n=0; n<nsp; n++){
      A.set(0.);
      for (int i=0; i<nsp; i++) // evaluate at each solution point
	for (int j=0; j<order+1; j++)
	  A(0,i) += pow(r[i],j)*lc[n*nsp+j];
      for (int i=0; i<nsp; i++){
	if (fabs(A(0,i)) < 5.e-13) A(0,i) = 0.;
	std::cout << A(0,i) << " ";
      }
      std::cout << std::endl;
    }}

  A.deallocate();
}
