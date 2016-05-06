#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::update(const int& step,
				const int& stage)
{
  /*
  // this is one way to update interior nodes
  double a=rka(stage),ap,au,b=0.,c,eps=1.e-14;
  if (step > 0) b = bdf(0);
  for(int n=0; n<nNode-nNodeBd; n++){
    ap = 1./(a*dt(n));
    au = b/max(dtUnsteady,eps);
    c  = 1./(ap+au);
    for (int k=0; k<nq; k++) q(n,k) = c*(ap*qn(n,k)+au*q(n,k)-r(n,k)/v(n));
  }
  */


  // compute update for interior nodes
  Array2D<double> dq(nNode,nq);
  double a=rka(stage),ap,au,b=0.,c,eps=1.e-14;
  if (step > 0) b = bdf(0);
  for(int n=0; n<nNode-nNodeBd; n++){
    ap = 1./(a*dt(n));
    au = b/max(dtUnsteady,eps);
    c  = 1./(v(n)*(ap+au));
    for (int k=0; k<nq; k++) dq(n,k) = c*r(n,k);
  }


  /*
  // compute update for boundary nodes
  int m=0;
  double rb[nq],L[nq*nq],A[nq*nq];
  for (int n=nNode-nNodeBd; n<nNode; n++){
    ap = 1./(a*dt(n));
    au = b/max(dtUnsteady,eps);
    c  = v(n)*(ap+au);
    sys->lhsBCVectorSelfJacobian(1,&nodeBd(m),&ln(m,0),&q(n,0),&qa(n,0),&A[0]);
    sys->rhsBCSelectionMatrix(1,&nodeBd(m),&ln(m,0),&q(n,0),&qa(n,0),&L[0]);
    for (int k=0; k<nq*nq; k++) A[k] += c*L[k];
    matinv(nq,&A[0]);
    matmul(nq,nq,1,&A[0],&r(n,0),&dq(n,0));
    //for (int k=0; k<nq; k++) dq(n,k) *= .25;
    m++;
  }
  */


  // compute update for boundary nodes
  int m=0;
  double rb[nq],L[nq*nq],Ab[nq*nq],A[nq*nq],Ak[nq*nq],Am[nq*nq],qk[nq],qm[nq];
  for (int n=nNode-nNodeBd; n<nNode; n++){
    ap = 1./(a*dt(n));
    au = b/max(dtUnsteady,eps);
    c  = v(n)*(ap+au);
    sys->lhsBCVectorSelfJacobian(1,&nodeBd(m),&ln(m,0),&q(n,0),&qa(n,0),&Ab[0]);
    sys->rhsBCSelectionMatrix(1,&nodeBd(m),&ln(m,0),&q(n,0),&qa(n,0),&L[0]);
    for (int k=0; k<nq*nq; k++) A[k]  = Ab[k]+c*L[k];
    for (int k=0; k<nq*nq; k++) Ak[k] = ap*v(n)*L[k]-A[k];
    for (int k=0; k<nq*nq; k++) Am[k] = au*v(n)*L[k]+Ab[k];
    matmul(nq,nq,1,&Ak[0],&qn(n,0),&qk[0]);
    matmul(nq,nq,1,&Am[0],&q (n,0),&qm[0]);
    for (int k=0; k<nq; k++) qk[k] =-qk[k]-qm[k]+r(n,k);
    matinv(nq,&A[0]);
    matmul(nq,nq,1,&A[0],&qk[0],&dq(n,0));
    //for (int k=0; k<nq; k++) dq(n,k) *= .25;
    m++;
  }


  // smooth updates
  Array2D<double> dq0(nNode,nq),dqn(nNode,nq);
  Array1D<double> g(nNode);
  int n1,n2;
  double dqk,gn,gr;
  for(int n=0; n<nNode; n++)
    for (int k=0; k<nq; k++) dq0(n,k) = dq(n,k);
  g.set(0.);
  for (int n=0; n<nEdge; n++){
    n1     = edge(n,0);
    n2     = edge(n,1);
    g(n1) += smooth;
    g(n2) += smooth;
  }

  for (int l=0; l<2; l++){
    dqn.set(0.);
  for (int n=0; n<nEdge; n++){
    n1     = edge(n,0);
    n2     = edge(n,1);
    for (int k=0; k<nq; k++){
      dqk        = dq(n2,k)-dq(n1,k);
      dqn(n1,k) += dqk;
      dqn(n2,k) -= dqk;
    }}
  int nn=nNode;
  if (sourceMMS == 1) nn = nNode-nNodeBd;
  for (int n=0; n<nn; n++){
    gn = g(n);
    gr = 1./(1.+gn);
    for (int k=0; k<nq; k++)
      dq(n,k) =(dq0(n,k)+gn*dq(n,k)+smooth*dqn(n,k))*gr;
  }}


  // apply update to all nodes
  for(int n=0; n<nNode; n++)
    for (int k=0; k<nq; k++) q(n,k) = qn(n,k)-dq(n,k);


  // set additional variables
  sys->stepQAdd(nNode,
		&q(0,0),
		&qa(0,0));

  dq.deallocate();
  dq0.deallocate();
  dqn.deallocate();
  g.deallocate();
}
