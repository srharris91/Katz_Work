#include "sgetrf.h"
#include "sgetri.h"
#include "matinv.h"


void matinv(int& n,
	    double* A)
{
  int info,ipiv[n];
  double work[n];
  sgetrf_(n,n,&A[0],n,&ipiv[0],info);
  sgetri_(n,&A[0],n,&ipiv[0],&work[0],n,info);
}
