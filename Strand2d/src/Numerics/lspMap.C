#include "StrandBlockSolver.h"
#include "ssytrf.h"
#include "ssytri.h"
#include "ssycon.h"
#include "sgesvd.h"


void StrandBlockSolver::lspMap()
{
  int mm,nn,jj,ii,j1,j2,jm,jmax,nmax,info,ldu,ldvt,rows,cols,lwork;
  double dsm,dx,dy,ds,w,r1,r2,rs,cond,rcond,xn,yn,ax,ay,xu,xl,bu,bl,xcn,ycn;
  double b[4];
  double* sv;
  double* work1;
  double* uu;
  double* vt;
  double* dr;
  double** lspT;
  char u='U',jobu='n',jobvt='a';
  cond = 0.;
  ldu  = 1;
  ldvt = 2;
  cols = 2;

  for (int n=0; n<nNodes-nGnodes; n++){
    mm    = ncsp(n);
    rows  = mm;
    ii    = cols;
    if (rows < ii) ii = rows;
    jj    = cols;
    if (rows > jj) jj = rows;
    lwork = 1;
    if (3*ii+jj > lwork) lwork = 3*ii+jj;
    if (5*ii    > lwork) lwork = 5*ii;
    uu    = new double[ldu*ldu];
    sv    = new double[cols];
    vt    = new double[ldvt*cols];
    work1 = new double[lwork];
    dr    = new double[rows*2];
    lspT  = new double*[nPstr+1];
    for (int j=0; j<nPstr+1; j++) lspT[j] = new double[rows];
    rcond = 0.;
    info  = 0;
    for (int m=0; m<ldu*ldu;   m++) uu   [m] = 0.;
    for (int m=0; m<cols;      m++) sv   [m] = 0.;
    for (int m=0; m<ldvt*cols; m++) vt   [m] = 0.;
    for (int m=0; m<lwork;     m++) work1[m] = 0.;
    double work [2] = {0.,0.};
    double work2[4] = {0.,0.,0.,0.};
    int    iwork[2] = {0,0};
    int    ipiv [2] = {0,0};

  for (int j=1; j<nPstr+1; j++){ //mid strand nodes

    // coordinates of the mid-strand location in question
    jm = j-1;
    xn = .5*(x(0,j,n)+x(0,jm,n));
    yn = .5*(x(1,j,n)+x(1,jm,n));

    // find data centroid
    xcn = 0.;
    ycn = 0.;
    for (int m=0; m<mm; m++){
      nn  = csp[n][m];
      xcn = xcn+xc(0,j,nn);
      ycn = ycn+xc(1,j,nn);
    }
    xcn = xcn/double(mm);
    ycn = ycn/double(mm);

    // find plane which most closely fits surrounding cell centers
    for (int m=0; m<mm; m++){
      nn       = csp[n][m];
      dr[m   ] = xc(0,j,nn)-xcn;
      dr[m+mm] = xc(1,j,nn)-ycn;
    }
    sgesvd_(jobu,jobvt,rows,cols,dr,rows,sv,uu,ldu,vt,ldvt,work1,lwork,info);
    if (info != 0){
      cout << "\n*** svd procedure failure in lspMap ***" << endl;
      exit(0);
    }
    ax = vt[0];
    ay = vt[2];
    ds = 1./sqrt(ax*ax+ay*ay);
    ax = ax*ds;
    ay = ay*ds;

    // compute 2d least squares problem with projected distances
    // largest distance in stencil
    dsm = 0.;
    for (int m=0; m<mm; m++){
      nn  = csp[n][m];
      dx  = xc(0,j,nn)-xn;
      dy  = xc(1,j,nn)-yn;
      ds  = dx*ax+dy*ay;
      ds  = ds*ds;
      if (ds > dsm) dsm = ds;
    }
    dsm = 1./sqrt(dsm);

    // form least squares matrix
    jj  = jm;
    for (int m=0; m<4; m++) b[m] = 0.;
    for (int m=0; m<mm; m++){
      nn   = csp[n][m];
      dx  = xc(0,j,nn)-xn;
      dy  = xc(1,j,nn)-yn;
      ds   =(dx*ax+dy*ay)*dsm;
      w    = 1./(ds*ds);
      b[0] = b[0]+w;
      b[2] = b[2]+w*ds;
      b[3] = b[3]+w*ds*ds;
    }

    // find max abs row sum for condition number computation
    r1 = fabs(b[0])+fabs(b[2]);
    r2 = fabs(b[1])+fabs(b[3]);
                 rs = r1;
    if (r2 > rs) rs = r2;

    // invert matrix and determine condition number
    ii = 2;
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
    for (int m=0; m<mm; m++){
      nn         = csp[n][m];
      dx         = xc(0,j,nn)-xn;
      dy         = xc(1,j,nn)-yn;
      ds         =(dx*ax+dy*ay)*dsm;
      w          = 1./(ds*ds);
      lspT[j][m] =(b[0]*w    +
                   b[2]*w*ds);
    }
  }

  // interpolate projected coefficients to the nodal positions along each strand
  for (int j=0; j<nPstr+1; j++){
    if      (j == 0    ){
      j1 = 1;
      j2 = 2;
    }
    else if (j == nPstr){
      j1 = nPstr-1;
      j2 = nPstr;
    }
    else{
      j1 = j;
      j2 = j+1;
    }
    xn = xStr(j);
    jm = j1-1;
    xl = .5*(xStr(jm)+xStr(j1));
    jm = j2-1;
    xu = .5*(xStr(jm)+xStr(j2));
    bl =(xu-xn)/(xu-xl);
    bu =(xn-xl)/(xu-xl);
    indlsp(0,j,n,ii);
    for (int m=0; m<mm; m++){
      lsp[ii  ][m] = bl*lspT[j1][m];
      lsp[ii+1][m] = bu*lspT[j2][m];
    }
  }

  delete [] sv;
  delete [] uu;
  delete [] vt;
  delete [] work1;
  delete [] dr;
  for (int j=0; j<nPstr+1; j++) delete [] lspT[j];
  delete [] lspT;
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



  /*
  // try using volume averaging on outer boundary nodes
  for (int n=0; n<nNodes-nGnodes; n++){
    mm = ncsp(n);
  for (int j=nPstr; j<nPstr+1; j++){
    if      (j == 0    ) jj = 1;
    else if (j == nPstr) jj = nPstr-1;
    else                 jj = j;

    w  = 0.;
    for (int k=0; k<2; k++){
      indlsp(k,j,n,ii);
    for (int m=0; m<mm; m++){
      nn = csp[n][m];
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
  */
}
