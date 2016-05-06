#include "Strand2dFCSPTurbSA.h"


void Strand2dFCSPTurbSA::rhsInvFluxX(const int& npts,
			             const double* q,
			             const double* qa,
			             double* f)
{
  int iq,iqa;
  for (int n=0; n<npts; n++){
    iq      = nq *n;
    iqa     = nqa*n;
    f[iq  ] = qa[iqa+1]* q[iq  ]         ;
    f[iq+1] = qa[iqa+1]* q[iq+1]+qa[iqa] ;
    f[iq+2] = qa[iqa+1]* q[iq+2]         ;
    f[iq+3] = qa[iqa+1]*(q[iq+3]+qa[iqa]);
    f[iq+4] = qa[iqa+1]* q[iq+4]         ;
  }
}
