#include "Tri2dFCSPLam.h"


void Tri2dFCSPLam::rhsInvFluxY(const int& npts,
			       const double* q,
			       const double* qa,
			       double* g)
{
  int iq,iqa;
  for (int n=0; n<npts; n++){
    iq       = nq *n;
    iqa      = nqa*n;
    g[iq  ]  = qa[iqa+2]* q[iq  ]         ;
    g[iq+1]  = qa[iqa+2]* q[iq+1]         ;
    g[iq+2]  = qa[iqa+2]* q[iq+2]+qa[iqa] ;
    g[iq+3]  = qa[iqa+2]*(q[iq+3]+qa[iqa]);
  }
}
