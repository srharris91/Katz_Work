#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::rhsInviscid()
{
  // compute inviscid fluxes at each node in the mesh
    Array3D<double>
      f(nSurfNode,nStrandNode,nq),
      g(nSurfNode,nStrandNode,nq);
    sys->rhsInvFluxX(nSurfNode*nStrandNode,&q(0,0,0),&qa(0,0,0),&f(0,0,0));
    sys->rhsInvFluxY(nSurfNode*nStrandNode,&q(0,0,0),&qa(0,0,0),&g(0,0,0));


  // surface edge inviscid flux contributions
  if (surfOrder == 1 || surfOrder == 2){ //linear scheme
    int i1,i2,n1,n2;
    double xn1,yn1,xn2,yn2,f1,f2,fI;
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<nElemEdge; i++){
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
	    f1         = yn1*f(n1,j,k)-xn1*g(n1,j,k);
	    f2         = yn2*f(n2,j,k)-xn2*g(n2,j,k);
	    fI         = .5*(f1+f2);
	    r(n1,j,k) += fI;
	    r(n2,j,k) -= fI;
	  }}}
    for (int n=0; n<nBndNode; n++){
      n1 = bndNode(n);
      for (int j=0; j<nStrandNode; j++){
	xn1 = xn(n1,j);
	yn1 = yn(n1,j);
	for (int k=0; k<nq; k++){
	  fI         = yn1*f(n1,j,k)-xn1*g(n1,j,k);
	  r(n1,j,k) += bndSign(n)*fI;
	}}}
  }

  else if (surfOrder == 3){ //corrected scheme
    Array4D<double>
      fx(nSurfNode,nStrandNode,nq,2),
      gx(nSurfNode,nStrandNode,nq,2);
    gradient(f,fx);
    gradient(g,gx);
    int i1,i2,n1,n2;
    double xn1,yn1,xn2,yn2,xs1,ys1,xs2,ys2,xns1,yns1,xns2,yns2,
      fs1,fs2,f1,f2,fI,ds5=.5*deltaS;
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<nElemEdge; i++){
	i1 = elemEdge(i,0);
	i2 = elemEdge(i,1);
	n1 = surfElem(n,i1);
	n2 = surfElem(n,i2);
	for (int j=0; j<nStrandNode; j++){
	  xn1  = xn(n1,j);
	  yn1  = yn(n1,j);
	  xn2  = xn(n2,j);
	  yn2  = yn(n2,j);
	  xs1  = xs(n,i1,j);
	  ys1  = ys(n,i1,j);
	  xs2  = xs(n,i2,j);
	  ys2  = ys(n,i2,j);
	  xns1 = xns(n,i1,j);
	  yns1 = yns(n,i1,j);
	  xns2 = xns(n,i2,j);
	  yns2 = yns(n,i2,j);
	  for (int k=0; k<nq; k++){
	    fs1        = yn1*(xs1*fx(n1,j,k,0)+ys1*fx(n1,j,k,1))
	                -xn1*(xs1*gx(n1,j,k,0)+ys1*gx(n1,j,k,1));
	    fs2        = yn2*(xs2*fx(n2,j,k,0)+ys2*fx(n2,j,k,1))
	                -xn2*(xs2*gx(n2,j,k,0)+ys2*gx(n2,j,k,1));
	    fs1       += yns1*f(n1,j,k)-xns1*g(n1,j,k);
	    fs2       += yns2*f(n2,j,k)-xns2*g(n2,j,k);
	    f1         = yn1*f(n1,j,k)-xn1*g(n1,j,k);
	    f2         = yn2*f(n2,j,k)-xn2*g(n2,j,k);
	    f1        += ds5*fs1;
	    f2        -= ds5*fs2;
	    fI         = .5*(f1+f2);
	    r(n1,j,k) += fI;
	    r(n2,j,k) -= fI;
	  }}}
    for (int n=0; n<nBndNode; n++){
      n1 = bndNode(n);
      for (int j=0; j<nStrandNode; j++){
	xn1 = xn(n1,j);
	yn1 = yn(n1,j);
	for (int k=0; k<nq; k++){
	  fI         = yn1*f(n1,j,k)-xn1*g(n1,j,k);
	  r(n1,j,k) += bndSign(n)*fI;
	}}}
    fx.deallocate();
    gx.deallocate();
  }


  // strand edge inviscid flux contributions
  Array2D<double> fhat(nStrandNode,nq);
  int ni,j,jH,mH,i1,i2,n1,n2;
  double dnr=1./deltaN,ds5=.5*deltaS,ds8=.125*deltaS*deltaS,ds6=deltaS/6.,
    df[nq],s1,s2,sa1,sa2;
  Array4D<double> rs(nSurfElem,meshOrder+1,nStrandNode,nq);
  for (int n=0; n<nSurfElem; n++){
    for (int i=0; i<meshOrder+1; i++){
      ni = surfElem(n,i);

      // compute directed fluxes at all nodes along this strand
      for (int j=0; j<nStrandNode; j++)
	for (int k=0; k<nq; k++)
	  fhat(j,k) =-ys(n,i,j)*f(ni,j,k)+xs(n,i,j)*g(ni,j,k);
    
      // flux derivatives at root boundary nodes
      for (int j=0; j<nIcbNode; j++){
	for (int k=0; k<nq; k++) df[k] = 0.;
	for (int m=0; m<nIcb; m++)
	  for (int k=0; k<nq; k++)
	    df[k] += icb(j,m)*fhat(m,k);
	for (int k=0; k<nq; k++) rs(n,i,j,k) = dnr*df[k];
      }

      // flux derivatives at interior nodes
      for (int j=nIcbNode; j<nStrandNode-nIcbNode; j++){
	for (int k=0; k<nq; k++) df[k] = 0.;
	for (int m=0; m<nIci; m++)
	  for (int k=0; k<nq; k++)
	    df[k] += ici(m)*(fhat(j+m+1,k)-fhat(j-m-1,k));
	for (int k=0; k<nq; k++) rs(n,i,j,k) = dnr*df[k];
      }

      // flux derivatives at tip boundary nodes
      jH = nIcbNode-1;
      for (int j=nStrandNode-nIcbNode; j<nStrandNode; j++){
	for (int k=0; k<nq; k++) df[k] = 0.;
	mH = nIcb-1;
	for (int m=0; m<nIcb; m++){
	  for (int k=0; k<nq; k++)
	    df[k] -= icb(jH,mH)*fhat(nStrandNode-nIcb+m,k);
	  mH--;
	}
	for (int k=0; k<nq; k++) rs(n,i,j,k) = dnr*df[k];
	jH--;
      }}}

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
	    r(n1,j,k) += .5*deltaS*rs(n,i1,j,k);
	    r(n2,j,k) += .5*deltaS*rs(n,i2,j,k);
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
	    s1         = rs(n,i1,j,k);
	    s2         = rs(n,i2,j,k);
	    sa1        = s1+s2;
	    sa2        = s1+s2;
	    r(n1,j,k) += ds6*(s1+sa1);
	    r(n2,j,k) += ds6*(s2+sa2);
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
	      ss(i,j,k) += ls(i,m)*rs(n,m,j,k);
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
	    s1         = rs(n,i1,j,k);
	    s2         = rs(n,i2,j,k);
	    sa1        = s1+s2-ds5*(ss (i1,j,k)+ss (i2,j,k))
	                      -ds8*(s2s(i1,j,k)+s2s(i2,j,k));
	    sa2        = s1+s2+ds5*(ss (i1,j,k)+ss (i2,j,k))
	                      -ds8*(s2s(i1,j,k)+s2s(i2,j,k));
	    r(n1,j,k) += ds6*(s1+sa1);
	    r(n2,j,k) += ds6*(s2+sa2);
	  }}}
    ss.deallocate();
    s2s.deallocate();
  }


  // clean up
  f.deallocate();
  g.deallocate();
  fhat.deallocate();
  rs.deallocate();
}
