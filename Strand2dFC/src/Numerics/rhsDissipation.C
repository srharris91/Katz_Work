#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::rhsDissipation(const int& j)
{
  // compute gradient of Q in the surface direction
  gradient(q,qx,j);


  // compute limiter
  limit(j);


  // surface edge dissipation contributions
  int i1,i2,n1,n2;
  double xs1,ys1,xn1,yn1,xs2,ys2,xn2,yn2,qs1,qs2,Ax,Ay,
    qL[nq],qR[nq],f[nq],dq[nq],ds5=.5*deltaS;
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<nElemEdge; i++){
      i1  = elemEdge(i,0);
      i2  = elemEdge(i,1);
      n1  = surfElem(n,i1);
      n2  = surfElem(n,i2);
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
	qs1   = qx(n1,j,k,0)*xs1+qx(n1,j,k,1)*ys1;
	qs2   = qx(n2,j,k,0)*xs2+qx(n2,j,k,1)*ys2;
	qL[k] = q(n1,j,k)+lim(n1,j,k)*ds5*qs1;
	qR[k] = q(n2,j,k)-lim(n2,j,k)*ds5*qs2;
	dq[k] = qR[k]-qL[k];
      }
      sys->rhsDisFlux(1,&Ax,&Ay,&qL[0],&qR[0],&dq[0],&f[0]);
      for (int k=0; k<nq; k++){
	d(n1,j,k) -= .5*f[k];
	d(n2,j,k) += .5*f[k];
      }}


  // strand edge dissipation contributions
  int ni,j1,nj;
  double dnr=.5/deltaN;
  Array3D<double> ds(nSurfElem,meshOrder+1,nq);
  ds.set(0.);
  if (j == 0) //root node
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<meshOrder+1; i++){
	ni = surfElem(n,i);
	j1 = dcn2[j][0];
	nj = dcn2[j][1];
	for (int k=0; k<nq; k++) dq[k] = 0.;
	for (int m=0; m<nj; m++){
	  for (int k=0; k<nq; k++) dq[k] += dcn1[j][m]*q(ni,j1,k);
	  j1++;
	}
	Ax =-.5*(ys(n,i,j)+ys(n,i,j+1));
	Ay = .5*(xs(n,i,j)+xs(n,i,j+1));
	sys->rhsDisFlux(1,&Ax,&Ay,&q(ni,j,0),&q(ni,j+1,0),&dq[0],&f[0]);
	for (int k=0; k<nq; k++) ds(n,i,k) += dnr*f[k];
	
	j1 = dcn2[j+1][0];
	nj = dcn2[j+1][1];
	for (int k=0; k<nq; k++) dq[k] = 0.;
	for (int m=0; m<nj; m++){
	  for (int k=0; k<nq; k++) dq[k] += dcn1[j+1][m]*q(ni,j1,k);
	  j1++;
	}
	Ax =-.5*(ys(n,i,j)+ys(n,i,j+1));
	Ay = .5*(xs(n,i,j)+xs(n,i,j+1));
	sys->rhsDisFlux(1,&Ax,&Ay,&q(ni,j,0),&q(ni,j+1,0),&dq[0],&f[0]);
	for (int k=0; k<nq; k++) ds(n,i,k) -= dnr*f[k];
      }
      
  else if (j == nStrandNode-1) //tip node
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<meshOrder+1; i++){
	ni = surfElem(n,i);
	j1 = dcn2[j][0];
	nj = dcn2[j][1];
	for (int k=0; k<nq; k++) dq[k] = 0.;
	for (int m=0; m<nj; m++){
	  for (int k=0; k<nq; k++) dq[k] += dcn1[j][m]*q(ni,j1,k);
	  j1++;
	}
	Ax =-.5*(ys(n,i,j-1)+ys(n,i,j));
	Ay = .5*(xs(n,i,j-1)+xs(n,i,j));
	sys->rhsDisFlux(1,&Ax,&Ay,&q(ni,j-1,0),&q(ni,j,0),&dq[0],&f[0]);
	for (int k=0; k<nq; k++) ds(n,i,k) += dnr*f[k];
	
	j1 = dcn2[j+1][0];
	nj = dcn2[j+1][1];
	for (int k=0; k<nq; k++) dq[k] = 0.;
	for (int m=0; m<nj; m++){
	  for (int k=0; k<nq; k++) dq[k] += dcn1[j+1][m]*q(ni,j1,k);
	  j1++;
	}
	Ax =-.5*(ys(n,i,j-1)+ys(n,i,j));
	Ay = .5*(xs(n,i,j-1)+xs(n,i,j));
	sys->rhsDisFlux(1,&Ax,&Ay,&q(ni,j-1,0),&q(ni,j,0),&dq[0],&f[0]);
	for (int k=0; k<nq; k++) ds(n,i,k) -= dnr*f[k];
      }
  
  else //interior nodes
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<meshOrder+1; i++){
	ni = surfElem(n,i);
	j1 = dcn2[j][0];
	nj = dcn2[j][1];
	for (int k=0; k<nq; k++) dq[k] = 0.;
	for (int m=0; m<nj; m++){
	  for (int k=0; k<nq; k++) dq[k] += dcn1[j][m]*q(ni,j1,k);
	  j1++;
	}
	Ax =-.5*(ys(n,i,j-1)+ys(n,i,j));
	Ay = .5*(xs(n,i,j-1)+xs(n,i,j));
	sys->rhsDisFlux(1,&Ax,&Ay,&q(ni,j-1,0),&q(ni,j,0),&dq[0],&f[0]);
	for (int k=0; k<nq; k++) ds(n,i,k) += dnr*f[k];
	
	j1 = dcn2[j+1][0];
	nj = dcn2[j+1][1];
	for (int k=0; k<nq; k++) dq[k] = 0.;
	for (int m=0; m<nj; m++){
	  for (int k=0; k<nq; k++) dq[k] += dcn1[j+1][m]*q(ni,j1,k);
	  j1++;
	}
	Ax =-.5*(ys(n,i,j)+ys(n,i,j+1));
	Ay = .5*(xs(n,i,j)+xs(n,i,j+1));
	sys->rhsDisFlux(1,&Ax,&Ay,&q(ni,j,0),&q(ni,j+1,0),&dq[0],&f[0]);
	for (int k=0; k<nq; k++) ds(n,i,k) -= dnr*f[k];
      }


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
  ds.deallocate();
}
