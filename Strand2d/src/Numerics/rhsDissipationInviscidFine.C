#include "StrandBlockSolver.h"


void StrandBlockSolver::rhsDissipationInviscidFine()
{
  int c1,c2,n1,n2,jm,jp,fc,npts=1;
  double xq,yq,dxl,dyl,dxr,dyr,ql[nq],qr[nq],fi[nq],fd[nq];


  // unstructured faces
  for (int n=0; n<nEdges; n++){
    c1            = edge(0,n);
    c2            = edge(1,n);
    n1            = edgn(n);
    fc            = max(fClip(c1),fClip(c2));
    for (int j=1; j<fc+1; j++){
      jm          = j-1;
      xq          = .5*(x(0,jm,n1)+x(0,j,n1));
      yq          = .5*(x(1,jm,n1)+x(1,j,n1));
      dxl         = xq-xc(0,j,c1);
      dyl         = yq-xc(1,j,c1);
      dxr         = xq-xc(0,j,c2);
      dyr         = yq-xc(1,j,c2);
      for (int k=0; k<nq; k++){
	ql[k]     = q(k,j,c1)+(dxl*qx(k,0,j,c1)+dyl*qx(k,1,j,c1))*lims(k,j,n);
	qr[k]     = q(k,j,c2)+(dxr*qx(k,0,j,c2)+dyr*qx(k,1,j,c2))*lims(k,j,n);
      }
      //ql[4] = q(4,j,c1);
      //qr[4] = q(4,j,c2);
      sys->rhsInvFlux(npts,&facs(0,j,n),&xvs(j,n),&ql[0],&qr[0],&fi[0]);
      sys->rhsDisFlux(npts,&facs(0,j,n),&xvs(j,n),&ql[0],&qr[0],&fd[0]);
      for (int k=0; k<nq; k++){
	r(k,j,c1) += (fi[k]-fd[k]);
	r(k,j,c2) -= (fi[k]-fd[k]);
      }}
  }


  // structured faces
  for (int n=0; n<nFaces-nGfaces; n++){
    n1             = face(0,n);
    n2             = face(1,n);
    for (int j=0; j<fClip(n)+1; j++){
      jp           = j+1;
      xq           = .5*(x(0,j,n1)+x(0,j,n2));
      yq           = .5*(x(1,j,n1)+x(1,j,n2));
      dxl          = xq-xc(0,j ,n);
      dyl          = yq-xc(1,j ,n);
      dxr          = xq-xc(0,jp,n);
      dyr          = yq-xc(1,jp,n);
      for (int k=0; k<nq; k++){
	ql[k]      = q(k,j ,n)+(dxl*qx(k,0,j ,n)+dyl*qx(k,1,j ,n))*limu(k,j,n);
	qr[k]      = q(k,jp,n)+(dxr*qx(k,0,jp,n)+dyr*qx(k,1,jp,n))*limu(k,j,n);
      }
      //ql[4] = q(4,j ,n);
      //qr[4] = q(4,jp,n);
      sys->rhsInvFlux(npts,&facu(0,j,n),&xvu(j,n),&ql[0],&qr[0],&fi[0]);
      sys->rhsDisFlux(npts,&facu(0,j,n),&xvu(j,n),&ql[0],&qr[0],&fd[0]);
      for (int k=0; k<nq; k++){
	r(k,j ,n) += (fi[k]-fd[k]);
	r(k,jp,n) -= (fi[k]-fd[k]);
      }}
  }
}
