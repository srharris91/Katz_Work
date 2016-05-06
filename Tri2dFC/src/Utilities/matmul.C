void matmul(const int& d1,
	    const int& d2,
	    const int& d3,
	    const double* B,
	    const double* C,
	    double* A)
{
  int m,n;
  for (int i=0; i<d1; i++)
  for (int j=0; j<d3; j++){
    m    = i*d3+j;
    n    = i*d2;
    A[m] = 0.;
    for (int k=0; k<d2; k++) A[m] += B[n+k]*C[k*d3+j];
  }
}
