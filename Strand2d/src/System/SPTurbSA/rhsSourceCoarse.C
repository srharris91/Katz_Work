#include "StrandSPTurbSA.h"
#include <math.h>

void StrandSPTurbSA::rhsSourceCoarse(const int& npts,
				     const double* v,
				     const double* q,
				     const double* qa,
				     double* r)
{
  int iq;
  for (int n=0; n<npts; n++){
    iq      = nq*n;
    r[iq  ] = 0.;
    r[iq+1] = 0.;
    r[iq+2] = 0.;
    r[iq+3] = 0.;
    r[iq+4] = 0.;
  }
}
