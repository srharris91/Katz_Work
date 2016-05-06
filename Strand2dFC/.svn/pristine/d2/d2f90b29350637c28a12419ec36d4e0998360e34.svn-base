#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::rhsViscous(const int& j)
{
  // gradients of qa in n-direction and s-direction
  int j1,nj,nm;
  double dnr=1./deltaN;
  Array2D<double> qan(nSurfNode,nqaGradQa);
  qan.set(0.);
  for (int n=0; n<nSurfNode; n++){
    j1 = icn2[j][0];
    nj = icn2[j][1];
    for (int m=0; m<nj; m++){
      for (int k=0; k<nqaGradQa; k++)
	qan(n,k) += dnr*icn1[j][m]*qa(n,j1,iqagrad(k));
      j1++;
    }}
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<meshOrder+1; i++){ // ith point in the element
      j1 = icn2[j][0];
      nj = icn2[j][1];
      for (int mm=0; mm<nj; mm++){
	for (int k=0; k<nqaGradQa; k++) qas(n,i,j1,k) = 0.;
	j1++;
      }
      for (int m=0; m<meshOrder+1; m++){ // mth Lagrange poly. in mapping
	nm = surfElem(n,m);
	j1 = icn2[j][0];
	nj = icn2[j][1];
	for (int mm=0; mm<nj; mm++){
	  for (int k=0; k<nqaGradQa; k++)
	    qas(n,i,j1,k) += ls(i,m)*qa(nm,j1,iqagrad(k));
	  j1++;
	}}}


  // surface edge viscous flux contributions
  int ni,i1,i2,n1,n2;
  double jac1,xs1,ys1,xn1,yn1,f1,f2,fI,
    qax[nqaGradQa],qay[nqaGradQa],fv[nq],gv[nq];
  Array2D<double> fhat(meshOrder+1,nq);
  if (surfOrder == 1 || surfOrder == 2) //linear scheme
    for (int n=0; n<nSurfElem; n++){
      for (int i=0; i<meshOrder+1; i++){ // transform to x-y derivatives
	                                 // and compute viscous fluxes
	ni   = surfElem(n,i);
	jac1 = 1./jac(n,i,j);
	xs1  = xs(n,i,j)*jac1;
	ys1  = ys(n,i,j)*jac1;
	xn1  = xn(ni,j)*jac1;
	yn1  = yn(ni,j)*jac1;
	for (int k=0; k<nqaGradQa; k++){
	  qax[k] = yn1*qas(n,i,j,k)-ys1*qan(ni,k);
	  qay[k] =-xn1*qas(n,i,j,k)+xs1*qan(ni,k);
	}
	sys->rhsVisFlux(1,&q(ni,j,0),&qa(ni,j,0),&qax[0],&qay[0],
			&fv[0],&gv[0]);
	xn1  = xn(ni,j);
	yn1  = yn(ni,j);
	for (int k=0; k<nq; k++) fhat(i,k) = yn1*fv[k]-xn1*gv[k];
      }
      for (int i=0; i<nElemEdge; i++){ // flux balance
	i1 = elemEdge(i,0);
	i2 = elemEdge(i,1);
	n1 = surfElem(n,i1);
	n2 = surfElem(n,i2);
	for (int k=0; k<nq; k++){
	  f1         = fhat(i1,k);
	  f2         = fhat(i2,k);
	  fI         = .5*(f1+f2);
	  d(n1,j,k) -= fI;
	  d(n2,j,k) += fI;
	}}}

  else if (surfOrder == 3){ //corrected scheme
    double ds5=.5*deltaS;
    Array2D<double> fs(meshOrder+1,nq);
    for (int n=0; n<nSurfElem; n++){
      for (int i=0; i<meshOrder+1; i++){ // transform to x-y derivatives
	                                 // and compute viscous fluxes
	ni   = surfElem(n,i);
	jac1 = 1./jac(n,i,j);
	xs1  = xs(n,i,j)*jac1;
	ys1  = ys(n,i,j)*jac1;
	xn1  = xn(ni,j)*jac1;
	yn1  = yn(ni,j)*jac1;
	for (int k=0; k<nqaGradQa; k++){
	  qax[k] = yn1*qas(n,i,j,k)-ys1*qan(ni,k);
	  qay[k] =-xn1*qas(n,i,j,k)+xs1*qan(ni,k);
	}
	sys->rhsVisFlux(1,&q(ni,j,0),&qa(ni,j,0),&qax[0],&qay[0],
			&fv[0],&gv[0]);
	xn1  = xn(ni,j);
	yn1  = yn(ni,j);
	for (int k=0; k<nq; k++) fhat(i,k) = yn1*fv[k]-xn1*gv[k];
      }

      fs.set(0.);
      for (int i=0; i<meshOrder+1; i++) // s-derivative of fhat
	for (int m=0; m<meshOrder+1; m++) // mth Lagrange poly. in mapping
	  for (int k=0; k<nq; k++)
	    fs(i,k) += ls(i,m)*fhat(m,k);

      for (int i=0; i<nElemEdge; i++){ // flux balance
	i1 = elemEdge(i,0);
	i2 = elemEdge(i,1);
	n1 = surfElem(n,i1);
	n2 = surfElem(n,i2);
	for (int k=0; k<nq; k++){
	  f1         = fhat(i1,k)+ds5*fs(i1,k);
	  f2         = fhat(i2,k)-ds5*fs(i2,k);
	  fI         = .5*(f1+f2);
	  d(n1,j,k) -= fI;
	  d(n2,j,k) += fI;
	}}}
    fs.deallocate();
  }

  // s-boundary node contributions to viscous flux and penalty
  int en,ei;
  double Pdnr=1.,uw[2],pbv[nq]; //for now, use fixed walls
  uw[0] = 0.;
  uw[1] = 0.;
  for (int n=0; n<nBndNode; n++){ //boundary nodes
    n1 = bndNode(n);
    en = bndElem(n,0);
    ei = bndElem(n,1);

    // transform to x-y derivatives and compute viscous fluxes
    jac1 = 1./jac(en,ei,j);
    xs1  = xs(en,ei,j)*jac1;
    ys1  = ys(en,ei,j)*jac1;
    xn1  = xn(n1,j)*jac1;
    yn1  = yn(n1,j)*jac1;
    for (int k=0; k<nqaGradQa; k++){
      qax[k] = yn1*qas(en,ei,j,k)-ys1*qan(n1,k);
      qay[k] =-xn1*qas(en,ei,j,k)+xs1*qan(n1,k);
    }
    sys->rhsVisFlux(1,&q(n1,j,0),&qa(n1,j,0),&qax[0],&qay[0],
		    &fv[0],&gv[0]);
    xn1  = xn(n1,j);
    yn1  = yn(n1,j);
    for (int k=0; k<nq; k++){ //viscous flux
      fv[k]      = yn1*fv[k]-xn1*gv[k];
      d(n1,j,k) -= bndSign(n)*fv[k];
    }

    for (int k=0; k<nq; k++){ //viscous penalty
      gv[k]      = yn1*bndDataVis(n,j,k,0)-xn1*bndDataVis(n,j,k,1);
      pbv[k]     = gv[k]-fv[k];
      d(n1,j,k) -= bndSign(n)*pbv[k];
    }}


  // strand edge viscous flux contributions
  Array3D<double> ds(nSurfElem,meshOrder+1,nq);
  Array2D<double> bn(nq,nqaGradQa),bnj(nq,nqaGradQa);
  ds.set(0.);
  int l1,nl;
  double df[nq],gshat[nq],gnhat[nq],qag[nqaGradQa],dnrr=dnr/deltaN;
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<meshOrder+1; i++){ // ith point in the element
      ni = surfElem(n,i);
      
      // s-portion of viscous flux derivatives along this strand
      for (int k=0; k<nq; k++) df[k] = 0.;
      j1 = icn2[j][0];
      nj = icn2[j][1];
      for (int m=0; m<nj; m++){
	jac1 = jac(n,i,j1);
	xs1  = xs(n,i,j1);
	ys1  = ys(n,i,j1);
	xn1  = xn(ni,j1);
	yn1  = yn(ni,j1);
	sys->rhsVisFluxS(1,&jac1,&xs1,&ys1,&xn1,&yn1,&q(ni,j1,0),&qa(ni,j1,0),
			 &qas(n,i,j1,0),&gshat[0]);
	for (int k=0; k<nq; k++) df[k] += icn1[j][m]*gshat[k];
	j1++;
      }
      for (int k=0; k<nq; k++) ds(n,i,k) =-dnr*df[k];

      // n-portion of viscous flux derivatives along this strand
      for (int k=0; k<nq; k++) df[k] = 0.;
      j1 = vcn2[j][0];
      nj = vcn2[j][1];
      for (int m=0; m<nj; m++){
	bnj.set(0.);
	l1 = vcn3[j][m][0];
	nl = vcn3[j][m][1];
	for (int l=0; l<nl; l++){
	  jac1 = jac(n,i,l1);
	  xs1  = xs(n,i,l1);
	  ys1  = ys(n,i,l1);
	  xn1  = xn(ni,l1);
	  yn1  = yn(ni,l1);
	  sys->rhsVisFluxNCoeff(1,&jac1,&xs1,&ys1,&xn1,&yn1,
				&q(ni,l1,0),&qa(ni,l1,0),&bn(0,0));
	  for (int ii=0; ii<nq; ii++)
	    for (int jj=0; jj<nqaGradQa; jj++)
	      bnj(ii,jj) += vcn1[j][m][l]*bn(ii,jj);
	  l1++;
	}
	for (int k=0; k<nqaGradQa; k++) qag[k] = qa(ni,j1,iqagrad(k));
	matmul(nq,nqaGradQa,1,&bnj(0,0),&qag[0],&gnhat[0]);
	for (int k=0; k<nq; k++) df[k] += gnhat[k];
	j1++;
      }
      for (int k=0; k<nq; k++) ds(n,i,k) -= dnrr*df[k];
    }


  // viscous flux penalties for strand roots
  Pdnr = Pinv0/deltaN;
  double qan1[nqaGradQa];
  if (j == 0)
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<meshOrder+1; i++){ // ith point in the element
	ni = surfElem(n,i);

	// s-viscous flux
	jac1 = jac(n,i,j);
	xs1  = xs(n,i,j);
	ys1  = ys(n,i,j);
	xn1  = xn(ni,j);
	yn1  = yn(ni,j);
	sys->rhsVisFluxS(1,&jac1,&xs1,&ys1,&xn1,&yn1,&q(ni,j,0),
			 &qa(ni,j,0),&qas(n,i,j,0),&gshat[0]);

	// n-viscous flux
	sys->rhsVisFluxNCoeff(1,&jac1,&xs1,&ys1,&xn1,&yn1,
			      &q(ni,j,0),&qa(ni,j,0),&bn(0,0));
	for (int k=0; k<nqaGradQa; k++) qan1[k] = 0.;
	for (int m=0; m<nBndVis; m++)
	  for (int k=0; k<nqaGradQa; k++)
	    qan1[k] += vcn4[m]*qa(ni,m,iqagrad(k));
	for (int k=0; k<nqaGradQa; k++) qan1[k] /= deltaN;
	matmul(nq,nqaGradQa,1,&bn(0,0),&qan1[0],&gnhat[0]);

	for (int k=0; k<nq; k++){
	  gv[k]      =-ys1*surfDataVis(ni,0,k,0)+xs1*surfDataVis(ni,0,k,1);
	  fv[k]      = gshat[k]+gnhat[k];
	  pbv[k]     = gv[k]-fv[k];
	  ds(n,i,k) += Pdnr*pbv[k];
	}}


  // viscous flux penalties for strand tips
  if (j == nStrandNode-1)
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<meshOrder+1; i++){ // ith point in the element
	ni = surfElem(n,i);

	// s-viscous flux
	jac1 = jac(n,i,j);
	xs1  = xs(n,i,j);
	ys1  = ys(n,i,j);
	xn1  = xn(ni,j);
	yn1  = yn(ni,j);
	sys->rhsVisFluxS(1,&jac1,&xs1,&ys1,&xn1,&yn1,&q(ni,j,0),
			 &qa(ni,j,0),&qas(n,i,j,0),&gshat[0]);

	// n-viscous flux
	sys->rhsVisFluxNCoeff(1,&jac1,&xs1,&ys1,&xn1,&yn1,
			      &q(ni,j,0),&qa(ni,j,0),&bn(0,0));
	for (int k=0; k<nqaGradQa; k++) qan1[k] = 0.;
	for (int m=0; m<nBndVis; m++)
	  for (int k=0; k<nqaGradQa; k++)
	    qan1[k] -= vcn4[m]*qa(ni,nStrandNode-1-m,iqagrad(k));
	for (int k=0; k<nqaGradQa; k++) qan1[k] /= deltaN;
	matmul(nq,nqaGradQa,1,&bn(0,0),&qan1[0],&gnhat[0]);

	for (int k=0; k<nq; k++){
	  gv[k]      =-ys1*surfDataVis(ni,1,k,0)+xs1*surfDataVis(ni,1,k,1);
	  fv[k]      = gshat[k]+gnhat[k];
	  pbv[k]     = gv[k]-fv[k];
	  ds(n,i,k) -= Pdnr*pbv[k];
	}}


  // source treatment
  int ne;
  double ns;
  for (int n=0; n<nSurfNode; n++)
    for (int i=psp2S(n); i<psp2S(n+1); i++){
      ne = psp1S(i,0);
      ni = psp1S(i,1);
      ns = wsp1S(i);
      for (int k=0; k<nq; k++) d(n,j,k) += ns*ds(ne,ni,k);
    }


  // clean up
  fhat.deallocate();
  qan.deallocate();
  ds.deallocate();
  bn.deallocate();
  bnj.deallocate();
}
