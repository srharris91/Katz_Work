#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::update(const int& step,
				   const int& stage,
				   const int& j)
{
  // compute update for all nodes
  if (stage == 0) pseudoTime(j);
  Array2D<double> A(nq,nq),dq(nSurfNode,nq);
  double a=rka(stage),b=0.,c,eps=1.e-14,dqm[nq],ap,au;
  if (step > 0) b = bdf(0);
  au = b/max(dtUnsteady,eps);
  for (int n=0; n<nSurfNode; n++){
    ap = 1./(a*dt(n,j));
    c  = v(n,j)*(ap+au);
    for (int k=0; k<nq; k++) dqm[k]  =-au*v(n,j)*(q(n,j,k)-qn(n,j,k));
    for (int k=0; k<nq; k++) dqm[k] += r(n,j,k);
    for (int k=0; k<nq; k++)
      for (int l=0; l<nq; l++) A(k,l) = Adl(n,k,l);
    for (int k=0; k<nq; k++) A(k,k) += c;
    matinv(nq,&A(0,0));
    matmul(nq,nq,1,&A(0,0),&dqm[0],&dq(n,0));
  }


  // compute update for strand root nodes (at this time, I treat viscous
  // wall strongly, which requires this older treatment)
  if (j == 0){
    double L[nq*nq],Aa[nq*nq],Ab[nq*nq],Am[nq*nq],dqk[nq];
    for (int n=0; n<nSurfNode; n++){
      ap = 1./(a*dt(n,j));
      c  = v(n,j)*(ap+au);
      sys->lhsBCVectorSelfJacobian(1,&surfNodeTag(n,0),&sn(n,0,0),
				   &q(n,j,0),&qa(n,j,0),&Ab[0]);
      sys->rhsBCSelectionMatrix(1,&surfNodeTag(n,0),&sn(n,0,0),
				&q(n,j,0),&qa(n,j,0),&L[0]);
      for (int k=0; k<nq; k++)
	for (int l=0; l<nq; l++) A(k,l) = Adl(n,k,l);
      for (int k=0; k<nq; k++) A(k,k) += c;
      matmul(nq,nq,nq,&L[0],&A(0,0),&Aa[0]);
      for (int k=0; k<nq*nq; k++) Aa[k] += Ab[k];
      for (int k=0; k<nq*nq; k++) Am[k]  =-(au*v(n,j)*L[k]+Ab[k]);
      for (int k=0; k<nq;    k++) dqm[k] = q(n,j,k)-qn(n,j,k);
      matmul(nq,nq,1,&Am[0],&dqm[0],&dqk[0]);
      for (int k=0; k<nq; k++) dqk[k] += r(n,j,k);
      matinv(nq,&Aa[0]);
      matmul(nq,nq,1,&Aa[0],&dqk[0],&dq(n,0));
    }}

  /*
  if (j == 0 || j == nStrandNode-1)
    for (int n=0; n<nSurfNode; n++)
      for (int k=0; k<nq; k++) dq(n,k) = 0.;
  */


  // smooth updates in the unstructured direction with Jacobi iterations
  Array2D<double> dq0(nSurfNode,nq),dqn(nSurfNode,nq);
  Array1D<double> g(nSurfNode);
  int i1,i2,n1,n2;
  double dqq,gn,gr,smoothSurf=.5*smooth;
  for (int n=0; n<nSurfNode; n++)
    for (int k=0; k<nq; k++) dq0(n,k) = dq(n,k);
  g.set(0.);
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<nElemEdge; i++){
      i1     = elemEdge(i,0);
      i2     = elemEdge(i,1);
      n1     = surfElem(n,i1);
      n2     = surfElem(n,i2);
      g(n1) += smoothSurf;
      g(n2) += smoothSurf;
  }
  for (int l=0; l<2; l++){
    dqn.set(0.);
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<nElemEdge; i++){
	i1 = elemEdge(i,0);
	i2 = elemEdge(i,1);
	n1 = surfElem(n,i1);
	n2 = surfElem(n,i2);
	for (int k=0; k<nq; k++){
	  dqq        = dq(n2,k)-dq(n1,k);
	  dqn(n1,k) += dqq;
	  dqn(n2,k) -= dqq;
	}}
    for (int n=0; n<nSurfNode; n++){
      gn = g(n);
      gr = 1./(1.+gn);
      for (int k=0; k<nq; k++)
	dq(n,k) =(dq0(n,k)+gn*dq(n,k)+smoothSurf*dqn(n,k))*gr;
    }}


  // apply update to all nodes
  for (int n=0; n<nSurfNode; n++)
    for (int k=0; k<nq; k++) q(n,j,k) = qn(n,j,k)-dq(n,k);


  // set additional variables
  for (int n=0; n<nSurfNode; n++) sys->stepQAdd(1,&q(n,j,0),&qa(n,j,0));


  // clean up
  A.deallocate();
  dq.deallocate();
  dq0.deallocate();
  dqn.deallocate();
  g.deallocate();
}
