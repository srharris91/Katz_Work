#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::lhsViscous(const int& j)
{
  // viscous variable Jacobian
  Array3D<double> G(nSurfNode,nqaGradQa,nq);
  for (int n=0; n<nSurfNode; n++)
    sys->lhsVisVariableJacobian(1,&q(n,j,0),&qa(n,j,0),&G(n,0,0));


  // strand edge viscous flux Jacobian contributions
  int ni,mj,l1,nl,j1,nj,lj;
  double jac1,xs1,ys1,xn1,yn1,dnrr=1./deltaN/deltaN,dqn[nq],dqan[nqa];
  Array2D<double> A(nq,nq),Aj(nq,nq),bn(nq,nqaGradQa),bnj(nq,nqaGradQa);
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<meshOrder+1; i++){
      ni = surfElem(n,i);

      // derivative with respect to viscous variables
      bn.set(0.);
      mj = vcn2[j][2];
      l1 = vcn3[j][mj][0];
      nl = vcn3[j][mj][1];
      for (int l=0; l<nl; l++){
	jac1 = jac(n,i,l1);
	xs1  = xs(n,i,l1);
	ys1  = ys(n,i,l1);
	xn1  = xn(ni,l1);
	yn1  = yn(ni,l1);
	sys->rhsVisFluxNCoeff(1,&jac1,&xs1,&ys1,&xn1,&yn1,
			      &q(ni,l1,0),&qa(ni,l1,0),&bnj(0,0));
	for (int ii=0; ii<nq; ii++)
	  for (int jj=0; jj<nqaGradQa; jj++)
	    bn(ii,jj) += vcn1[j][mj][l]*bnj(ii,jj);
	l1++;
      }
      matmul(nq,nqaGradQa,nq,&bn(0,0),&G(ni,0,0),&A(0,0));

      // derivative with respect to coefficient matrix
      for (int k=0; k<nq; k++) dqn[k] = 0.;
      for (int k=0; k<nqa; k++) dqan[k] = 0.;
      j1 = vcn2[j][0];
      nj = vcn2[j][1];
      for (int m=0; m<nj; m++){
	lj = vcn3[j][m][2];
	for (int k=0; k<nq; k++) dqn[k] += vcn1[j][m][lj]*q(ni,j1,k);
	for (int k=0; k<nqa; k++) dqan[k] += vcn1[j][m][lj]*qa(ni,j1,k);
	j1++;
      }
      jac1 = jac(n,i,j);
      xs1  = xs(n,i,j);
      ys1  = ys(n,i,j);
      xn1  = xn(ni,j);
      yn1  = yn(ni,j);
      sys->lhsVisFluxNCoeffJacobian(1,&jac1,&xs1,&ys1,&xn1,&yn1,
				    &q(ni,j,0),&qa(ni,j,0),
				    &dqn[0],&dqan[0],&Aj(0,0));
      for (int k=0; k<nq; k++)
	for (int l=0; l<nq; l++) A(k,l) += Aj(k,l);

      // add to LHS Jacobian
      for (int k=0; k<nq; k++)
	for (int l=0; l<nq; l++) Ads(n,i,k,l) -= dnrr*A(k,l);
    }


  // strand root penalty terms
  dnrr *= Pinv0;
  if (j == 0)
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<meshOrder+1; i++){
	ni = surfElem(n,i);
	
	// derivative with respect to viscous variables
	jac1 = jac(n,i,j);
	xs1  = xs(n,i,j);
	ys1  = ys(n,i,j);
	xn1  = xn(ni,j);
	yn1  = yn(ni,j);
	sys->rhsVisFluxNCoeff(1,&jac1,&xs1,&ys1,&xn1,&yn1,
			      &q(ni,j,0),&qa(ni,j,0),&bn(0,0));
	matmul(nq,nqaGradQa,nq,&bn(0,0),&G(ni,0,0),&A(0,0));
	for (int k=0; k<nq; k++)
	  for (int l=0; l<nq; l++) A(k,l) *= vcn4[0];
	
	// derivative with respect to coefficient matrix
	for (int k=0; k<nq; k++) dqn[k] = 0.;
	for (int k=0; k<nqa; k++) dqan[k] = 0.;
	for (int m=0; m<nBndVis; m++){
	  for (int k=0; k<nq; k++)
	    dqn[k] += vcn4[m]*q(ni,m,k);
	  for (int k=0; k<nqa; k++)
	    dqan[k] += vcn4[m]*qa(ni,m,k);
	}
	sys->lhsVisFluxNCoeffJacobian(1,&jac1,&xs1,&ys1,&xn1,&yn1,
				      &q(ni,j,0),&qa(ni,j,0),
				      &dqn[0],&dqan[0],&Aj(0,0));
	for (int k=0; k<nq; k++)
	  for (int l=0; l<nq; l++) A(k,l) += Aj(k,l);
	
	// add to LHS Jacobian
	for (int k=0; k<nq; k++)
	  for (int l=0; l<nq; l++) Ads(n,i,k,l) -= dnrr*A(k,l);
      }


  // strand tip penalty terms
  if (j == nStrandNode-1)
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<meshOrder+1; i++){
	ni = surfElem(n,i);
	
	// derivative with respect to viscous variables
	jac1 = jac(n,i,j);
	xs1  = xs(n,i,j);
	ys1  = ys(n,i,j);
	xn1  = xn(ni,j);
	yn1  = yn(ni,j);
	sys->rhsVisFluxNCoeff(1,&jac1,&xs1,&ys1,&xn1,&yn1,
			      &q(ni,j,0),&qa(ni,j,0),&bn(0,0));
	matmul(nq,nqaGradQa,nq,&bn(0,0),&G(ni,0,0),&A(0,0));
	for (int k=0; k<nq; k++)
	  for (int l=0; l<nq; l++) A(k,l) *=(-vcn4[0]);
	
	// derivative with respect to coefficient matrix
	for (int k=0; k<nq; k++) dqn[k] = 0.;
	for (int k=0; k<nqa; k++) dqan[k] = 0.;
	for (int m=0; m<nBndVis; m++){
	  for (int k=0; k<nq; k++)
	    dqn[k] -= vcn4[m]*q(ni,nStrandNode-1-m,k);
	  for (int k=0; k<nqa; k++)
	    dqan[k] -= vcn4[m]*qa(ni,nStrandNode-1-m,k);
	}
	sys->lhsVisFluxNCoeffJacobian(1,&jac1,&xs1,&ys1,&xn1,&yn1,
				      &q(ni,j,0),&qa(ni,j,0),
				      &dqn[0],&dqan[0],&Aj(0,0));
	for (int k=0; k<nq; k++)
	  for (int l=0; l<nq; l++) A(k,l) += Aj(k,l);
	
	// add to LHS Jacobian
	for (int k=0; k<nq; k++)
	  for (int l=0; l<nq; l++) Ads(n,i,k,l) += dnrr*A(k,l);
      }


  // clean up
  G.deallocate();
  A.deallocate();
  Aj.deallocate();
  bn.deallocate();
  bnj.deallocate();
}
