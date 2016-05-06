#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::rhsViscous()
{
  // gradients of qa in n-direction
  int jH,mH;
  double dnr=1./deltaN;
  Array3D<double>
    fv(meshOrder+1,nStrandNode,nq),
    gv(meshOrder+1,nStrandNode,nq),
    qas(meshOrder+1,nStrandNode,nqaGradQa);
  Array3D<double> qan(nSurfNode,nStrandNode,nqaGradQa);
  qan.set(0.);
  for (int n=0; n<nSurfNode; n++){
    for (int j=0; j<nIcbNode; j++){ //root boundary nodes
      for (int m=0; m<nIcb; m++)
	for (int k=0; k<nqaGradQa; k++)
	  qan(n,j,k) += icb(j,m)*qa(n,m,iqagrad(k));
      for (int k=0; k<nqaGradQa; k++) qan(n,j,k) *= dnr;
    }
    for (int j=nIcbNode; j<nStrandNode-nIcbNode; j++){ //interior nodes
      for (int m=0; m<nIci; m++)
	for (int k=0; k<nqaGradQa; k++)
	  qan(n,j,k) += ici(m)*(qa(n,j+m+1,iqagrad(k))-qa(n,j-m-1,iqagrad(k)));
      for (int k=0; k<nqaGradQa; k++) qan(n,j,k) *= dnr;
    }
    jH = nIcbNode-1;
    for (int j=nStrandNode-nIcbNode; j<nStrandNode; j++){ //tip boundary nodes
      mH = nIcb-1;
      for (int m=0; m<nIcb; m++){
	for (int k=0; k<nqaGradQa; k++)
	  qan(n,j,k) -= icb(jH,mH)*qa(n,nStrandNode-nIcb+m,iqagrad(k));
	mH--;
      }
      for (int k=0; k<nqaGradQa; k++) qan(n,j,k) *= dnr;
      jH--;
    }}


  // surface edge viscous flux contributions
  if (surfOrder == 1 || surfOrder == 2){ //linear scheme
    int nm,ni,i1,i2,n1,n2,en,ei;
    double jac1,xs1,ys1,xn1,yn1,xn2,yn2,f1,f2,fI,qax[nqaGradQa],qay[nqaGradQa];
    for (int n=0; n<nSurfElem; n++){

      qas.set(0.); // local gradients of qa in s-direction
      for (int i=0; i<meshOrder+1; i++) // ith point in the element
	for (int m=0; m<meshOrder+1; m++){ // mth Lagrange poly. in mapping
	  nm = surfElem(n,m);
	  for (int j=0; j<nStrandNode; j++)
	    for (int k=0; k<nqaGradQa; k++)
	      qas(i,j,k) += ls(i,m)*qa(nm,j,iqagrad(k));
	}

      for (int i=0; i<meshOrder+1; i++){ // transform to x-y derivatives
	                                 // and compute viscous fluxes
	ni = surfElem(n,i);
	for (int j=0; j<nStrandNode; j++){
	  jac1 = 1./jac(n,i,j);
	  xs1  = xs(n,i,j)*jac1;
	  ys1  = ys(n,i,j)*jac1;
	  xn1  = xn(ni,j)*jac1;
	  yn1  = yn(ni,j)*jac1;
	  for (int k=0; k<nqaGradQa; k++){
	    qax[k] = yn1*qas(i,j,k)-ys1*qan(ni,j,k);
	    qay[k] =-xn1*qas(i,j,k)+xs1*qan(ni,j,k);
	  }
	  sys->rhsVisFlux(1,&q(ni,j,0),&qa(ni,j,0),&qax[0],&qay[0],
			  &fv(i,j,0),&gv(i,j,0));
	}}

      for (int i=0; i<nElemEdge; i++){ // flux balance
	i1 = elemEdge(i,0);
	i2 = elemEdge(i,1);
	n1 = surfElem(n,i1);
	n2 = surfElem(n,i2);
	for (int j=0; j<nStrandNode; j++){
	  xn1 = xn(n1,j);
	  yn1 = yn(n1,j);
	  xn2 = xn(n2,j);
	  yn2 = yn(n2,j);
	  for (int k=0; k<nq; k++){
	    f1         = yn1*fv(i1,j,k)-xn1*gv(i1,j,k);
	    f2         = yn2*fv(i2,j,k)-xn2*gv(i2,j,k);
	    fI         = .5*(f1+f2);
	    d(n1,j,k) -= fI;
	    d(n2,j,k) += fI;
	  }}}}

    for (int n=0; n<nBndNode; n++){ //boundary nodes
      n1 = bndNode(n);
      en = bndElem(n,0);
      ei = bndElem(n,1);

      qas.set(0.); // local gradients of qa in s-direction
      for (int m=0; m<meshOrder+1; m++){ // mth Lagrange poly. in mapping
	nm = surfElem(en,m);
	for (int j=0; j<nStrandNode; j++)
	  for (int k=0; k<nqaGradQa; k++)
	    qas(ei,j,k) += ls(ei,m)*qa(nm,j,iqagrad(k));
      }

      // transform to x-y derivatives and compute viscous fluxes
      for (int j=0; j<nStrandNode; j++){
	jac1 = 1./jac(en,ei,j);
	xs1  = xs(en,ei,j)*jac1;
	ys1  = ys(en,ei,j)*jac1;
	xn1  = xn(n1,j)*jac1;
	yn1  = yn(n1,j)*jac1;
	for (int k=0; k<nqaGradQa; k++){
	  qax[k] = yn1*qas(ei,j,k)-ys1*qan(n1,j,k);
	  qay[k] =-xn1*qas(ei,j,k)+xs1*qan(n1,j,k);
	}
	sys->rhsVisFlux(1,&q(n1,j,0),&qa(n1,j,0),&qax[0],&qay[0],
			&fv(ei,j,0),&gv(ei,j,0));
      }

      for (int j=0; j<nStrandNode; j++){
	xn1 = xn(n1,j);
	yn1 = yn(n1,j);
	for (int k=0; k<nq; k++){
	  fI         = yn1*fv(ei,j,k)-xn1*gv(ei,j,k);
	  d(n1,j,k) -= bndSign(n)*fI;
	}}}
  }

  else if (surfOrder == 3){ //corrected scheme
    cout << "***\nrhsViscous not yet done for surfOrder = 3.***" << endl;
    exit(0);
  }


  // strand edge viscous flux contributions
  Array4D<double> ds(nSurfElem,meshOrder+1,nStrandNode,nq);
  Array3D<double> bn(nStrandNode,nq,nqaGradQa);
  Array2D<double> qaS(nStrandNode,nqaGradQa),gshat(nStrandNode,nq),
    qaGradQa(nStrandNode,nqaGradQa),bnj(nq,nqaGradQa);
  ds.set(0.);
  int nm,ni,i1,i2,n1,n2,m0=(nVci-1)/2;
  double ds5=.5*deltaS,ds8=.125*deltaS*deltaS,ds6=deltaS/6.,
    jac1,xs1,ys1,xn1,yn1,s1,s2,sa1,sa2,df[nq],gnhat[nq];
  for (int n=0; n<nSurfElem; n++){
    //for (int n=0; n<0; n++){

    for (int i=0; i<meshOrder+1; i++){ // ith point in the element
      ni = surfElem(n,i);
      qaS.set(0.); // local gradients of qa in s-direction
      for (int m=0; m<meshOrder+1; m++){ // mth Lagrange poly. in mapping
	nm = surfElem(n,m);
	for (int j=0; j<nStrandNode; j++)
	  for (int k=0; k<nqaGradQa; k++)
	    qaS(j,k) += ls(i,m)*qa(nm,j,iqagrad(k));
      }

      // compute s-portion of viscous fluxes at all nodes along this strand
      for (int j=0; j<nStrandNode; j++){
	jac1 = jac(n,i,j);
	xs1  = xs(n,i,j);
	ys1  = ys(n,i,j);
	xn1  = xn(ni,j);
	yn1  = yn(ni,j);
	sys->rhsVisFluxS(1,&jac1,&xs1,&ys1,&xn1,&yn1,&q(ni,j,0),&qa(ni,j,0),
			 &qaS(j,0),&gshat(j,0));
      }

      // s-portion of viscous flux derivatives at root boundary nodes
      for (int j=0; j<nIcbNode; j++){
	for (int k=0; k<nq; k++) df[k] = 0.;
	for (int m=0; m<nIcb; m++)
	  for (int k=0; k<nq; k++)
	    df[k] += icb(j,m)*gshat(m,k);
	for (int k=0; k<nq; k++) ds(n,i,j,k) =-dnr*df[k];
      }

      // s-portion of viscous flux derivatives at interior nodes
      for (int j=nIcbNode; j<nStrandNode-nIcbNode; j++){
	for (int k=0; k<nq; k++) df[k] = 0.;
	for (int m=0; m<nIci; m++)
	  for (int k=0; k<nq; k++)
	    df[k] += ici(m)*(gshat(j+m+1,k)-gshat(j-m-1,k));
	for (int k=0; k<nq; k++) ds(n,i,j,k) =-dnr*df[k];
      }

      // s-portion of viscous flux derivatives at tip boundary nodes
      jH = nIcbNode-1;
      for (int j=nStrandNode-nIcbNode; j<nStrandNode; j++){
	for (int k=0; k<nq; k++) df[k] = 0.;
	mH = nIcb-1;
	for (int m=0; m<nIcb; m++){
	  for (int k=0; k<nq; k++)
	    df[k] -= icb(jH,mH)*gshat(nStrandNode-nIcb+m,k);
	  mH--;
	}
	for (int k=0; k<nq; k++) ds(n,i,j,k) =-dnr*df[k];
	jH--;
      }

      // compute coefficient matrix of n-portion of viscous fluxes
      // at all nodes along this strand, along with relevant qa variables
      for (int j=0; j<nStrandNode; j++){
	jac1 = jac(n,i,j);
	xs1  = xs(n,i,j);
	ys1  = ys(n,i,j);
	xn1  = xn(ni,j);
	yn1  = yn(ni,j);
	sys->rhsVisFluxNCoeff(1,&jac1,&xs1,&ys1,&xn1,&yn1,
			      &q(ni,j,0),&qa(ni,j,0),&bn(j,0,0));
	for (int k=0; k<nqaGradQa; k++) qaGradQa(j,k) = qa(ni,j,iqagrad(k));
      }

      // n-portion of viscous flux derivatives at root boundary nodes
      for (int j=0; j<nVcbNode; j++){

	for (int k=0; k<nq; k++) df[k] = 0.;

	for (int m=0; m<nVcb2; m++){
	  for (int ii=0; ii<nq; ii++)
	    for (int jj=0; jj<nqaGradQa; jj++) bnj(ii,jj) = 0.;
	  for (int mm=vcbIndex(j,m,0); mm<=vcbIndex(j,m,1); mm++)
	    for (int ii=0; ii<nq; ii++)
	      for (int jj=0; jj<nqaGradQa; jj++)
		bnj(ii,jj) += vcb(j,m,mm-vcbIndex(j,m,0))*bn(j+mm,ii,jj);

	  matmul(nq,nqaGradQa,1,&bnj(0,0),&qaGradQa(m,0),&gnhat[0]);
	  for (int k=0; k<nq; k++) df[k] += gnhat[k];
	}
	for (int k=0; k<nq; k++) ds(n,i,j,k) -= dnr*dnr*df[k];
      }

      //  n-portion of viscous flux derivatives at interior nodes
      for (int j=nVcbNode; j<nStrandNode-nVcbNode; j++){

	for (int k=0; k<nq; k++) df[k] = 0.;

	for (int m=0; m<nVci; m++){
	  for (int ii=0; ii<nq; ii++)
	    for (int jj=0; jj<nqaGradQa; jj++) bnj(ii,jj) = 0.;
	  for (int mm=vciIndex(m,0); mm<=vciIndex(m,1); mm++)
	    for (int ii=0; ii<nq; ii++)
	      for (int jj=0; jj<nqaGradQa; jj++)
		bnj(ii,jj) += vci(m,mm-vciIndex(m,0))*bn(j+mm,ii,jj);

	  matmul(nq,nqaGradQa,1,&bnj(0,0),&qaGradQa(j-m0+m,0),&gnhat[0]);
	  for (int k=0; k<nq; k++) df[k] += gnhat[k];
	}
	for (int k=0; k<nq; k++) ds(n,i,j,k) -= dnr*dnr*df[k];
      }

      //  n-portion of viscous flux derivatives at tip boundary nodes
      jH = nVcbNode-1;
      for (int j=nStrandNode-nVcbNode; j<nStrandNode; j++){

	for (int k=0; k<nq; k++) df[k] = 0.;
	mH = nVcb2-1;

	for (int m=0; m<nVcb2; m++){
	  for (int ii=0; ii<nq; ii++)
	    for (int jj=0; jj<nqaGradQa; jj++) bnj(ii,jj) = 0.;
	  for (int mm=vcbIndex(jH,mH,0); mm<=vcbIndex(jH,mH,1); mm++)
	    for (int ii=0; ii<nq; ii++)
	      for (int jj=0; jj<nqaGradQa; jj++)
		bnj(ii,jj) += vcb(jH,mH,mm-vcbIndex(jH,mH,0))
		              *bn(nStrandNode-nVcbNode-mm,ii,jj);

	  matmul(nq,nqaGradQa,1,&bnj(0,0),
		 &qaGradQa(nStrandNode-nVcb2+m,0),&gnhat[0]);
	  for (int k=0; k<nq; k++) df[k] += gnhat[k];
	  mH--;
	}
	for (int k=0; k<nq; k++) ds(n,i,j,k) -= dnr*dnr*df[k];
	jH--;
      }
    }
  }


  // source treatment
  if (surfOrder == 1 || surfOrder == 2) //mass lumped
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<nElemEdge; i++){
	i1 = elemEdge(i,0);
	i2 = elemEdge(i,1);
	n1 = surfElem(n,i1);
	n2 = surfElem(n,i2);
	for (int j=0; j<nStrandNode; j++)
	  for (int k=0; k<nq; k++){
	    d(n1,j,k) += .5*deltaS*ds(n,i1,j,k);
	    d(n2,j,k) += .5*deltaS*ds(n,i2,j,k);
	  }}

  /*
  else if (surfOrder == 2) //Galerkin
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<nElemEdge; i++){
	i1 = elemEdge(i,0);
	i2 = elemEdge(i,1);
	n1 = surfElem(n,i1);
	n2 = surfElem(n,i2);
	for (int j=0; j<nStrandNode; j++)
	  for (int k=0; k<nq; k++){
	    s1         = ds(n,i1,j,k);
	    s2         = ds(n,i2,j,k);
	    sa1        = s1+s2;
	    sa2        = s1+s2;
	    d(n1,j,k) += ds6*(s1+sa1);
	    d(n2,j,k) += ds6*(s2+sa2);
	  }}
  */

  else if (surfOrder == 3){ //corrected source by edge
    Array3D<double> ss(meshOrder+1,nStrandNode,nq);
    Array3D<double> s2s(meshOrder+1,nStrandNode,nq);
    for (int n=0; n<nSurfElem; n++){
      ss.set(0.);
      s2s.set(0.);
      for (int i=0; i<meshOrder+1; i++) // ith point in the element
	for (int m=0; m<meshOrder+1; m++) // mth Lagrange poly. in mapping
	  for (int j=0; j<nStrandNode; j++)
	    for (int k=0; k<nq; k++)
	      ss(i,j,k) += ls(i,m)*ds(n,m,j,k);
      for (int i=0; i<meshOrder+1; i++) // ith point in the element
	for (int m=0; m<meshOrder+1; m++) // mth Lagrange poly. in mapping
	  for (int j=0; j<nStrandNode; j++)
	    for (int k=0; k<nq; k++)
	      s2s(i,j,k) += ls(i,m)*ss(m,j,k);
      for (int i=0; i<nElemEdge; i++){
	i1 = elemEdge(i,0);
	i2 = elemEdge(i,1);
	n1 = surfElem(n,i1);
	n2 = surfElem(n,i2);
	for (int j=0; j<nStrandNode; j++)
	  for (int k=0; k<nq; k++){
	    s1         = ds(n,i1,j,k);
	    s2         = ds(n,i2,j,k);
	    sa1        = s1+s2-ds5*(ss (i1,j,k)+ss (i2,j,k))
	                      -ds8*(s2s(i1,j,k)+s2s(i2,j,k));
	    sa2        = s1+s2+ds5*(ss (i1,j,k)+ss (i2,j,k))
	                      -ds8*(s2s(i1,j,k)+s2s(i2,j,k));
	    d(n1,j,k) += ds6*(s1+sa1);
	    d(n2,j,k) += ds6*(s2+sa2);
	  }}}
    ss.deallocate();
    s2s.deallocate();
  }


  // clean up
  fv.deallocate();
  gv.deallocate();
  qas.deallocate();
  qan.deallocate();
  ds.deallocate();
  qaS.deallocate();
  gshat.deallocate();
  bn.deallocate();
  qaGradQa.deallocate();
  bnj.deallocate();
}
