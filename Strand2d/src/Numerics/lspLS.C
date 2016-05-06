#include "StrandBlockSolver.h"
#include "ssytrf.h"
#include "ssytri.h"
#include "ssycon.h"


void StrandBlockSolver::lspLS()
{
  int mm,nn,jj,ii,jm,jmax,nmax,info;
  double dsm,xci,yci,dx,dy,ds,w,r1,r2,r3,rs,cond,rcond,xn,yn;
  double b[9];
  char u='U';
  cond = 0.;
  rcond = 0.;
  info = 0;
  double work [3] = {0.,0.,0.};
  double work2[3] = {0.,0.,0.};
  int    iwork[3] = {0,0,0};
  int    ipiv [3] = {0,0,0};


  for (int n=0; n<nNodes-nGnodes; n++){
    mm = ncsp(n);
  for (int j=0; j<nPstr+1; j++){
    if      (j == 0    ) jm = 1;
    else if (j == nPstr) jm = nPstr-1;
    else                 jm = j;

    // coordinates of the central node in question
    xn = x(0,j,n);
    yn = x(1,j,n);

    // largest distance in stencil
    jj  = jm;
    dsm = 0.;
    for (int k=0; k<2; k++){
    for (int m=0; m<mm; m++){
      nn  = csp[n][m];
      xci = xc(0,jj,nn);
      yci = xc(1,jj,nn);
      dx  = xci-xn;
      dy  = yci-yn;
      ds  = dx*dx+dy*dy;
      if (ds > dsm) dsm = ds;
    }
    jj++;
    }
    dsm = 1./sqrt(dsm);

    // form least squares matrix
    jj  = jm;
    for (int m=0; m<9; m++) b[m] = 0.;
    for (int k=0; k<2; k++){
    for (int m=0; m<mm; m++){
      nn   = csp[n][m];
      xci  = xc(0,jj,nn);
      yci  = xc(1,jj,nn);
      dx   =(xci-xn)*dsm;
      dy   =(yci-yn)*dsm;
      w    = 1./(dx*dx+dy*dy);
      b[0] = b[0]+w;
      b[3] = b[3]+w*dx;
      b[6] = b[6]+w*dy;
      b[4] = b[4]+w*dx*dx;
      b[7] = b[7]+w*dx*dy;
      b[8] = b[8]+w*dy*dy;
    }
    jj++;
    }

    // find max abs row sum for condition number computation
    r1 = fabs(b[0])+fabs(b[3])+fabs(b[6]);
    r2 = fabs(b[3])+fabs(b[4])+fabs(b[7]);
    r3 = fabs(b[6])+fabs(b[7])+fabs(b[8]);
                 rs = r1;
    if (r2 > rs) rs = r2;
    if (r3 > rs) rs = r3;

    // invert matrix and determine condition number
    ii = 3;
    ssytrf_(u,ii,b,ii,ipiv,work,ii,info);
    ssycon_(u,ii,b,ii,ipiv,rs,rcond,work2,iwork,info);
    ssytri_(u,ii,b,ii,ipiv,work,info);
    rcond = 1./rcond;
    if (rcond > cond){
      nmax = n;
      jmax = j;
      cond = rcond;
    }
    
    // form and store least squares coefficient
    jj  = jm;
    for (int k=0; k<2; k++){
      indlsp(k,j,n,ii);
    for (int m=0; m<mm; m++){
      nn   = csp[n][m];
      xci  = xc(0,jj,nn);
      yci  = xc(1,jj,nn);
      dx   =(xci-xn)*dsm;
      dy   =(yci-yn)*dsm;
      w    = 1./(dx*dx+dy*dy);
      lsp[ii][m] =(b[0]*w    +
                   b[3]*w*dx +
                   b[6]*w*dy);
    }
    jj++;
    }
  }
  }

  // output condition information
  xn = x(0,jmax,nmax);
  yn = x(1,jmax,nmax);
  cout << "\nMaximum condition number for LS procedure: "
       << cond << endl
       << "Index of maximum condition number: "
       << nmax << "\t" << jmax << endl
       << "Coordinates of maximum condition number: "
       << xn << "\t" << yn << "\n" << endl;
}
