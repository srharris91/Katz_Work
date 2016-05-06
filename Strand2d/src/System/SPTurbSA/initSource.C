#include "StrandSPTurbSA.h"


void StrandSPTurbSA::initSource(const int& npts,
				const double* x,
				double* s)
{
  if (isolution == 3){
    cout << "\ninitSource not done yet in StrandSPTurbSA." << endl;
    exit(0);
  }
  else
    for (int n=0; n<nq*npts; n++) s[n] = 0.;
}


void StrandSPTurbSA::initSource(const int& npts,
				const int* tag,
				const double* x,
				double* s)
{
  if (isolution == 3){ //mms for boundary conditions not done yet
    cout << "\ninitSource not done yet in StrandSPTurbSA." << endl;
    exit(0);
  }
  else
    for (int n=0; n<nq*npts; n++) s[n] = 0.;
}
