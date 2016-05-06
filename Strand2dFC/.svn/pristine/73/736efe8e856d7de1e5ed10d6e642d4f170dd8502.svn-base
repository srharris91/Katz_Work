#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::lhsSource(const int& j)
{
  // gradients of qa in n-direction and s-direction
  int j1,nj,nm;
  double dnr=1./deltaN;
  Array2D<double> qan(nSurfNode,nqa);
  Array3D<double> qas(nSurfElem,meshOrder+1,nqa);
  qan.set(0.);
  for (int n=0; n<nSurfNode; n++){
    j1 = icn2[j][0];
    nj = icn2[j][1];
    for (int m=0; m<nj; m++){
      for (int k=0; k<nqa; k++)
	qan(n,k) += dnr*icn1[j][m]*qa(n,j1,k);
      j1++;
    }}
  qas.set(0.);
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<meshOrder+1; i++) // ith point in the element
      for (int m=0; m<meshOrder+1; m++){ // mth Lagrange poly. in mapping
	nm = surfElem(n,m);
	for (int k=0; k<nqa; k++)
	  qas(n,i,k) += ls(i,m)*qa(nm,j,k);
      }


  // source Jacobian computation
  int ni;
  double jac1,xs1,ys1,xn1,yn1,qax[nqa],qay[nqa];
  Array2D<double> A(nq,nq);
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<meshOrder+1; i++){
      ni   = surfElem(n,i);
      jac1 = 1./jac(n,i,j);
      xs1  = xs(n,i,j)*jac1;
      ys1  = ys(n,i,j)*jac1;
      xn1  = xn(ni,j)*jac1;
      yn1  = yn(ni,j)*jac1;
      for (int k=0; k<nqa; k++){
	qax[k] = yn1*qas(n,i,k)-ys1*qan(ni,k);
	qay[k] =-xn1*qas(n,i,k)+xs1*qan(ni,k);
      }
      sys->lhsSourceJacobian(1,&q(ni,j,0),&qa(ni,j,0),&qax[0],&qay[0],&A(0,0));
      for (int k=0; k<nq; k++)
	for (int l=0; l<nq; l++) Ads(n,i,k,l) -= A(k,l)/jac1;
    }


  // clean up
  qan.deallocate();
  qas.deallocate();
  A.deallocate();
}
