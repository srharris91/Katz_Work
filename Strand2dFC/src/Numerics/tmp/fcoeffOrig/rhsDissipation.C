#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::rhsDissipation()
{
  // compute gradient of Q in the surface direction
  gradient(q,qx);


  // compute limiter
  limit();


  // surface edge dissipation contributions
  int i1,i2,n1,n2;
  double xs1,ys1,xn1,yn1,xs2,ys2,xn2,yn2,qs1,qs2,Ax,Ay,
    qL[nq],qR[nq],f[nq],dqs[nq],ds5=.5*deltaS;
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<nElemEdge; i++){
      i1 = elemEdge(i,0);
      i2 = elemEdge(i,1);
      n1 = surfElem(n,i1);
      n2 = surfElem(n,i2);
      for (int j=0; j<nStrandNode; j++){
	xs1 = xs(n,i1,j);
	ys1 = ys(n,i1,j);
	xs2 = xs(n,i2,j);
	ys2 = ys(n,i2,j);
	xn1 = xn(n1,j);
	yn1 = yn(n1,j);
	xn2 = xn(n2,j);
	yn2 = yn(n2,j);
	Ax  = .5*(yn1+yn2);
	Ay  =-.5*(xn1+xn2);
	for (int k=0; k<nq; k++){
	  qs1    = qx(n1,j,k,0)*xs1+qx(n1,j,k,1)*ys1;
	  qs2    = qx(n2,j,k,0)*xs2+qx(n2,j,k,1)*ys2;
	  qL[k]  = q(n1,j,k)+lim(n1,j,k)*ds5*qs1;
	  qR[k]  = q(n2,j,k)-lim(n2,j,k)*ds5*qs2;
	  dqs[k] = qR[k]-qL[k];
	}
	sys->rhsDisFlux(1,&Ax,&Ay,&qL[0],&qR[0],&dqs[0],&f[0]);
	for (int k=0; k<nq; k++){
	  d(n1,j,k) -= .5*f[k];
	  d(n2,j,k) += .5*f[k];
	}}}


  // strand edge dissipation contributions
  int ni,jH,mH,m1=-nDci1/2,mh=-nDci/2;
  double dn5=.5*deltaN,dnr=.5/deltaN,ds8=.125*deltaS*deltaS,ds6=deltaS/6.,
    jac1,jac2,s1,s2,sa1,sa2,sL,sR;
  Array2D<double> dq1(nStrandNode+1,nq),dqh(nStrandNode+1,nq);
  Array4D<double> ds(nSurfElem,meshOrder+1,nStrandNode,nq);
  ds.set(0.);
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<meshOrder+1; i++){
      ni = surfElem(n,i);

      // store first order edge differences for this strand
      dq1.set(0.);
      for (int j=0; j<nDcbEdge1; j++) //root edges
	for (int m=0; m<nDcb1; m++)
	  for (int k=0; k<nq; k++)
	    dq1(j,k) += dcb1(j,m)*q(ni,m,k);
      for (int j=nDcbEdge1; j<nStrandNode-nDcbEdge1+1; j++) //interior edges
	for (int m=0; m<nDci1; m++)
	  for (int k=0; k<nq; k++) dq1(j,k) += dci1(m)*q(ni,j+m+m1,k);
      jH = nDcbEdge1-1;
      for (int j=nStrandNode-nDcbEdge1+1; j<nStrandNode+1; j++){ //tip edges
	mH = nDcb1-1;
	for (int m=0; m<nDcb1; m++){
	  for (int k=0; k<nq; k++)
	    dq1(j,k) -= dcb1(jH,mH)*q(ni,nStrandNode-nDcb1+m,k);
	  mH--;
	}
	jH--;
      }
      
      // store higher order edge differences for this strand
      dqh.set(0.);
      for (int j=0; j<nDcbEdge; j++) //root edges
	for (int m=0; m<nDcb; m++)
	  for (int k=0; k<nq; k++)
	    dqh(j,k) += dcb(j,m)*q(ni,m,k);
      for (int j=nDcbEdge; j<nStrandNode-nDcbEdge+1; j++) //interior edges
	for (int m=0; m<nDci; m++)
	  for (int k=0; k<nq; k++) dqh(j,k) += dci(m)*q(ni,j+m+mh,k);
      jH = nDcbEdge-1;
      for (int j=nStrandNode-nDcbEdge+1; j<nStrandNode+1; j++){ //tip edges
	mH = nDcb-1;
	for (int m=0; m<nDcb; m++){
	  for (int k=0; k<nq; k++)
	    dqh(j,k) -= dcb(jH,mH)*q(ni,nStrandNode-nDcb+m,k);
	  mH--;
	}
	jH--;
      }
      
      // compute limited edge differences for this strand
      for (int j=0; j<nStrandNode+1; j++)
	for (int k=0; k<nq; k++) dq1(j,k) += strandLim(ni,j,k)*dqh(j,k);
      
      int j=0; //root edge
      Ax =-.5*(ys(n,i,j)+ys(n,i,j+1));
      Ay = .5*(xs(n,i,j)+xs(n,i,j+1));
      sys->rhsDisFlux(1,&Ax,&Ay,&q(ni,j,0),&q(ni,j+1,0),&dq1(j,0),&f[0]);
      for (int k=0; k<nq; k++) ds(n,i,j  ,k) -= dnr*f[k];
      
      for (int j=1; j<nStrandNode; j++){ //interior edges
	Ax =-.5*(ys(n,i,j-1)+ys(n,i,j));
	Ay = .5*(xs(n,i,j-1)+xs(n,i,j));
	sys->rhsDisFlux(1,&Ax,&Ay,&q(ni,j-1,0),&q(ni,j,0),&dq1(j,0),&f[0]);
	for (int k=0; k<nq; k++){
	  ds(n,i,j-1,k) -= dnr*f[k];
	  ds(n,i,j  ,k) += dnr*f[k];
	}}
      
      j  = nStrandNode; //tip edge
      Ax =-.5*(ys(n,i,j-2)+ys(n,i,j-1));
      Ay = .5*(xs(n,i,j-2)+xs(n,i,j-1));
      sys->rhsDisFlux(1,&Ax,&Ay,&q(ni,j-2,0),&q(ni,j-1,0),&dq1(j,0),&f[0]);
      for (int k=0; k<nq; k++) ds(n,i,j-1,k) += dnr*f[k];
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
  dq1.deallocate();
  dqh.deallocate();
  ds.deallocate();
}
