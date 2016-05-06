#include <math.h>
#include "nchoosek.h"
#include "factorial.h"


double jacobiPoly(const int& n,
		  const int& alpha,
		  const int& beta,
		  const double& x)
{
  // See Canuto p. 70
  double P=0.;
  for (int i=0; i<=n; i++)
    P += nchoosek(n+alpha,i)*nchoosek(n+beta,n-i)*pow(x+1.,i)*pow(x-1.,n-i);
  P /= pow(2.,n);

  // ortho-normalize
  int i=alpha+beta;
  double a;
  a = pow(2.,i+1)/(double)(2*n+i+1)
    *(double)factorial(n+alpha)*(double)factorial(n+beta)
    /(double)factorial(n+i)/(double)factorial(n);
  P /= sqrt(a);

  return(P);
}
