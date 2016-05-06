#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::rhsInviscid(const int& j)
{
  // compute inviscid fluxes at each node in the mesh
  Array2D<double> f(nSurfNode,nq),g(nSurfNode,nq);
  for (int n=0; n<nSurfNode; n++){
    sys->rhsInvFluxX(1,&q(n,j,0),&qa(n,j,0),&f(n,0));
    sys->rhsInvFluxY(1,&q(n,j,0),&qa(n,j,0),&g(n,0));
  }


  // surface edge inviscid flux contributions
  int i1,i2,n1,n2;
  double xn1,yn1,xn2,yn2,f1,f2,fI;
  if (surfOrder == 1 || surfOrder == 2) //linear scheme
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<nElemEdge; i++){
	i1  = elemEdge(i,0);
	i2  = elemEdge(i,1);
	n1  = surfElem(n,i1);
	n2  = surfElem(n,i2);
	xn1 = xn(n1,j);
	yn1 = yn(n1,j);
	xn2 = xn(n2,j);
	yn2 = yn(n2,j);
	for (int k=0; k<nq; k++){
	  f1         = yn1*f(n1,k)-xn1*g(n1,k);
	  f2         = yn2*f(n2,k)-xn2*g(n2,k);
	  fI         = .5*(f1+f2);
	  r(n1,j,k) += fI;
	  r(n2,j,k) -= fI;
	}}

  else if (surfOrder == 3){ //corrected scheme
    Array3D<double> fx(nSurfNode,nq,2),gx(nSurfNode,nq,2);
    double xs1,ys1,xs2,ys2,xns1,yns1,xns2,yns2,fs1,fs2,ds5=.5*deltaS;

    int ni;
    double cx,cy;
    for (int n=0; n<nSurfNode; n++){
      for (int k=0; k<nq; k++){
	fx(n,k,0) = 0.;
	fx(n,k,1) = 0.;
	gx(n,k,0) = 0.;
	gx(n,k,1) = 0.;
      }
      for (int i=psp2(n); i<psp2(n+1); i++){
	ni = psp1(i);
	cx = gxc(i,j,0);
	cy = gxc(i,j,1);
	for (int k=0; k<nq; k++){
	  fx(n,k,0) += cx*f(ni,k);
	  fx(n,k,1) += cy*f(ni,k);
	  gx(n,k,0) += cx*g(ni,k);
	  gx(n,k,1) += cy*g(ni,k);
	}}}

    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<nElemEdge; i++){
	i1   = elemEdge(i,0);
	i2   = elemEdge(i,1);
	n1   = surfElem(n,i1);
	n2   = surfElem(n,i2);
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
	  fs1        = yn1*(xs1*fx(n1,k,0)+ys1*fx(n1,k,1))
	              -xn1*(xs1*gx(n1,k,0)+ys1*gx(n1,k,1));
	  fs2        = yn2*(xs2*fx(n2,k,0)+ys2*fx(n2,k,1))
	              -xn2*(xs2*gx(n2,k,0)+ys2*gx(n2,k,1));
	  fs1       += yns1*f(n1,k)-xns1*g(n1,k);
	  fs2       += yns2*f(n2,k)-xns2*g(n2,k);
	  f1         = yn1*f(n1,k)-xn1*g(n1,k);
	  f2         = yn2*f(n2,k)-xn2*g(n2,k);
	  f1        += ds5*fs1;
	  f2        -= ds5*fs2;
	  fI         = .5*(f1+f2);
	  r(n1,j,k) += fI;
	  r(n2,j,k) -= fI;
	}}
    fx.deallocate();
    gx.deallocate();
  }

  // s-boundary node contributions to boundary flux and penalty
  double uw[2],A[2],pb[nq],Pdnr=1.;
  uw[0] = 0.; //for now, use fixed walls in boundary conditions
  uw[1] = 0.;
  for (int n=0; n<nBndNode; n++){
    n1  = bndNode(n);
    xn1 = xn(n1,j);
    yn1 = yn(n1,j);
    for (int k=0; k<nq; k++){
      fI         = yn1*f(n1,k)-xn1*g(n1,k);
      r(n1,j,k) += bndSign(n)*fI;
    }

    A[0] = yn(n1,j);
    A[1] =-xn(n1,j);
    sys->rhsBCPenalty(1,&bndNodeTag(n),-(int)bndSign(n),&A[0],&Pdnr,
		      &q(n1,j,0),&qa(n1,j,0),&bndData(n,j,0),
		      &uw[0],&pb[0]);
    for (int k=0; k<nq; k++) r(n1,j,k) += bndSign(n)*pb[k];
  } 
  f.deallocate();
  g.deallocate();


  // strand edge inviscid flux contributions
  int ni,j1,nj;
  double dnr=1./deltaN,df[nq],fhat[nq],fs[nq],gs[nq];
  Array3D<double> rs(nSurfElem,meshOrder+1,nq);
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<meshOrder+1; i++){
      ni = surfElem(n,i);
      j1 = icn2[j][0];
      nj = icn2[j][1];
      for (int k=0; k<nq; k++) df[k] = 0.;
      for (int m=0; m<nj; m++){
	sys->rhsInvFluxX(1,&q(ni,j1,0),&qa(ni,j1,0),&fs[0]);
	sys->rhsInvFluxY(1,&q(ni,j1,0),&qa(ni,j1,0),&gs[0]);
	for (int k=0; k<nq; k++)
	  fhat[k] =-ys(n,i,j1)*fs[k]+xs(n,i,j1)*gs[k];
	for (int k=0; k<nq; k++) df[k] += icn1[j][m]*fhat[k];
	j1++;
      }
      for (int k=0; k<nq; k++) rs(n,i,k) = dnr*df[k];
    }

  if (j == 0){ // strand root boundary penalty
    Pdnr = Pinv0*dnr;
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<meshOrder+1; i++){
	ni   = surfElem(n,i);
	A[0] =-ys(n,i,j);
	A[1] = xs(n,i,j);
	sys->rhsBCPenalty(1,&surfNodeTag(ni,0),1,&A[0],&Pdnr,
			  &q(ni,j,0),&qa(ni,j,0),&surfData(ni,0,0),
			  &uw[0],&pb[0]);
	for (int k=0; k<nq; k++) rs(n,i,k) -= Pdnr*pb[k];
      }}
  
  if (j == nStrandNode-1){ // strand tip boundary penalty
    Pdnr = Pinv0*dnr;
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<meshOrder+1; i++){
	ni   = surfElem(n,i);
	A[0] =-ys(n,i,j);
	A[1] = xs(n,i,j);
	sys->rhsBCPenalty(1,&surfNodeTag(ni,1),-1,&A[0],&Pdnr,
			  &q(ni,j,0),&qa(ni,j,0),&surfData(ni,1,0),
			  &uw[0],&pb[0]);
	for (int k=0; k<nq; k++) rs(n,i,k) += Pdnr*pb[k];
      }}


  // source treatment
  int ne;
  double ns;
  for (int n=0; n<nSurfNode; n++)
    for (int i=psp2S(n); i<psp2S(n+1); i++){
      ne = psp1S(i,0);
      ni = psp1S(i,1);
      ns = wsp1S(i);
      for (int k=0; k<nq; k++) r(n,j,k) += ns*rs(ne,ni,k);
    }


  // clean up
  rs.deallocate();
}
