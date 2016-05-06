#include "StrandBlockSolver.h"


void StrandBlockSolver::lspVol()
{
  int mm,nn,jj,ii;
  double w;
  for (int n=0; n<nNodes-nGnodes; n++){
    mm = ncsp(n);
  for (int j=0; j<nPstr+1; j++){
    if      (j == 0    ) jj = 1;
    else if (j == nPstr) jj = nPstr-1;
    else                 jj = j;

    w  = 0.;
    for (int k=0; k<2; k++){
      indlsp(k,j,n,ii);
    for (int m=0; m<mm; m++){
      nn         = csp[n][m];
      lsp[ii][m] = v(jj,nn);
      w         += v(jj,nn);
    }
    jj++;
    }

    w = 1./w;
    for (int k=0; k<2; k++){
      indlsp(k,j,n,ii);
    for (int m=0; m<mm; m++){
      lsp[ii][m] = lsp[ii][m]*w;
    }}}}
}
