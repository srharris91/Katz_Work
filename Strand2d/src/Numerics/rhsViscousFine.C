#include "StrandBlockSolver.h"


void StrandBlockSolver::rhsViscousFine()
{
  int c1,c2,n1,n2,fc,jm,jp,npts=1;
  double dx1,dy1,dx2,dy2,ds,dq1,dq2,eps=1.e-14,
    qxe[nq],qye[nq],qaxe[nqa],qaye[nqa],qe[nq],qae[nqa],fv[nq];


  // unstructured faces
  for (int n=0; n<nEdges; n++){
    c1             = edge(0,n);
    c2             = edge(1,n);
    n1             = edgn(n);
    fc             = fClip(c1);
    if (fClip(c2) > fc) fc = fClip(c2);
    for (int j=1; j<fc+1; j++){
      jm           = j-1;
      dx1          = x (0,j,n1)-x (0,jm,n1);
      dy1          = x (1,j,n1)-x (1,jm,n1);
      dx2          = xc(0,j,c2)-xc(0,j ,c1);
      dy2          = xc(1,j,c2)-xc(1,j ,c1);
      ds           = 1./(dx1*dy2-dx2*dy1);
      for (int k=0; k<nq; k++){
	dq1        = qp(k,j,n1)-qp(k,jm,n1);
	dq2        = q (k,j,c2)-q (k,j ,c1);
	qxe[k]     = ds*( dy2*dq1-dy1*dq2);
	qye[k]     = ds*(-dx2*dq1+dx1*dq2);
	qe[k]      = .5*(qp (k,j,n1)+qp (k,jm,n1));
      }
      for (int k=0; k<nqa; k++){
	dq1        = qap(k,j,n1)-qap(k,jm,n1);
	dq2        = qa (k,j,c2)-qa (k,j ,c1);
	qaxe[k]    = ds*( dy2*dq1-dy1*dq2);
	qaye[k]    = ds*(-dx2*dq1+dx1*dq2);
	qae[k]     = .5*(qap(k,j,n1)+qap(k,jm,n1));
      }
      sys->rhsVisFlux(npts,&facs(0,j,n),&qe[0],&qae[0],&qxe[0],&qye[0],
		      &qaxe[0],&qaye[0],&fv[0]);
      for (int k=0; k<nq; k++){
	r(k,j,c1) -= fv[k];
	r(k,j,c2) += fv[k];
      }}}


  // structured faces
  for (int n=0; n<nFaces-nGfaces; n++){
    n1             = face(0,n);
    n2             = face(1,n);
    for (int j=0; j<fClip(n)+1; j++){
      jp           = j+1;
      dx1          = x (0,j ,n2)-x (0,j,n1);
      dy1          = x (1,j ,n2)-x (1,j,n1);
      dx2          = xc(0,jp,n )-xc(0,j,n );
      dy2          = xc(1,jp,n )-xc(1,j,n );
      ds           = dx1*dy2-dx2*dy1;
      if (fabs(ds) < eps) for (int k=0; k<nq; k++) fv[k] = 0.;
      else{
	ds         = 1./ds;
	for (int k=0; k<nq; k++){
	  dq1      = qp(k,j ,n2)-qp(k,j,n1);
	  dq2      = q (k,jp,n )-q (k,j,n );
	  qxe[k]   = ds*( dy2*dq1-dy1*dq2);
	  qye[k]   = ds*(-dx2*dq1+dx1*dq2);
	  qe[k]    = .5*(qp (k,j,n1)+qp (k,j,n2));
	}
	for (int k=0; k<nqa; k++){
	  dq1      = qap(k,j ,n2)-qap(k,j,n1);
	  dq2      = qa (k,jp,n )-qa (k,j,n );
	  qaxe[k]  = ds*( dy2*dq1-dy1*dq2);
	  qaye[k]  = ds*(-dx2*dq1+dx1*dq2);
	  qae[k]   = .5*(qap(k,j,n1)+qap(k,j,n2));
	}
	sys->rhsVisFlux(npts,&facu(0,j,n),&qe[0],&qae[0],&qxe[0],&qye[0],
			&qaxe[0],&qaye[0],&fv[0]);
      }
      for (int k=0; k<nq; k++){
	r(k,j ,n) -= fv[k];
	r(k,jp,n) += fv[k];
      }}}
}
