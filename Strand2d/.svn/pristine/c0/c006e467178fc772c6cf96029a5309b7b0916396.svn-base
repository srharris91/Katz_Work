#include "sgetrf.h"
#include "sgetri.h"
#include "matmul.h"


void nbtrluC(const int& il,
	     const int& ih,
	     const int& i,
	     const int& nq,
	     const double* a,
	     const double* b,
	     const double* c,
	     double* d,
	     double* e)
{
  int info,ipiv[nq];
  int n   = nq;
  int id  = nq*   (i  -il);
  int idm = nq*   (i-1-il);
  int ie  = nq*nq*(i  -il);
  int iem = nq*nq*(i-1-il);
  double aa[nq*nq],dd[nq],work[nq];
  if (i == il){
    for (int k=0; k<nq*nq; k++) aa[k] = b[k];
    sgetrf_(n,n,&aa[0],n,&ipiv[0],info);
    sgetri_(n,&aa[0],n,&ipiv[0],&work[0],n,info);
    matmul(n,n,1,&aa[0],&d[id],&dd[0]);
    for (int k=0; k<nq; k++) d[id+k] = dd[k];
    matmul(n,n,n,&aa[0],&c[0],&e[ie]);
  }
  else{
    matmul(n,n,n,&a[0],&e[iem],&aa[0]);
    for (int k=0; k<nq*nq; k++) aa[k] = b[k]-aa[k];
    sgetrf_(n,n,&aa[0],n,&ipiv[0],info);
    sgetri_(n,&aa[0],n,&ipiv[0],&work[0],n,info);
    matmul(n,n,1,&a[0],&d[idm],&dd[0]);
    for (int k=0; k<nq; k++) dd[k] = d[id+k]-dd[k];
    matmul(n,n,1,&aa[0],&dd[0],&d[id]);
    matmul(n,n,n,&aa[0],&c[0],&e[ie]);
  }
}


void nbtrbkC(const int& il,
	     const int& ih,
	     const int& nq,
	     double* d,
	     const double* e)
{
  int id,idp,ie,n=nq;
  double dd[nq];
  for (int i=ih-1; i>=il; i--){
    id  = nq*   (i  -il);
    idp = nq*   (i+1-il);
    ie  = nq*nq*(i  -il);
    matmul(n,n,1,&e[ie],&d[idp],&dd[0]);
    for (int k=0; k<nq; k++) d[id+k] -= dd[k];
  }
}
