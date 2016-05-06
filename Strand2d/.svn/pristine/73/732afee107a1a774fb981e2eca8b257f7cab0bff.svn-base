#include "StrandBlockSolver.h"


void StrandBlockSolver::gradQ(const int& mglevel)
{
  if (mglevel == 0 && gradient != 0 && gradQFlag == 0){

    // compute nodal values
    nodalQ(mglevel);


    // initialize gradients
    for (int n=0; n<nFaces+nBedges; n++)
      for (int j=0; j<nPstr+2; j++)
	for (int k=0; k<ndim; k++)
	  for (int m=0; m<nq; m++) qx(m,k,j,n) = 0.;


    // loop through unstructured edges
    int c1,c2,n1,n2,jm,jp,k;
    double Ax,Ay,sq;
    for (int n=0; n<nEdges; n++){
      c1            = edge(0,n);
      c2            = edge(1,n);
      n1            = edgn(n);
    for (int j=1; j<nPstr+1; j++){
      jm            = j-1;
      Ax            = facs(0,j,n);
      Ay            = facs(1,j,n);
    for (int kk=0; kk<nqGradQ; kk++){
      k             = iqgrad(kk);
      sq            = qp(k,jm,n1)+qp(k,j,n1);
      qx(k,0,j,c1) += Ax*sq;
      qx(k,1,j,c1) += Ay*sq;
      qx(k,0,j,c2) -= Ax*sq;
      qx(k,1,j,c2) -= Ay*sq;
    }}}


    // loop through structured edges
    for (int n=0; n<nFaces-nGfaces; n++){
      n1            = face(0,n);
      n2            = face(1,n);
    for (int j=0; j<nPstr+1; j++){
      jp            = j+1;
      Ax            = facu(0,j,n);
      Ay            = facu(1,j,n);
    for (int kk=0; kk<nqGradQ; kk++){
      k             = iqgrad(kk);
      sq            = qp(k,j,n1)+qp(k,j,n2);
      qx(k,0,j ,n) += Ax*sq;
      qx(k,1,j ,n) += Ay*sq;
      qx(k,0,jp,n) -= Ax*sq;
      qx(k,1,jp,n) -= Ay*sq;
    }}}


    // divide by twice the volume for the interior cells
    for (int n=0; n<nFaces-nGfaces; n++)
      for (int j=1; j<nPstr+1; j++)
	for (int m=0; m<ndim; m++)
	  for (int kk=0; kk<nqGradQ; kk++){
	    k            = iqgrad(kk);
	    qx(k,m,j,n) /= (2.*v(j,n));
	  }


    // surface, end, and boundary gradients
    double dx1,dy1,dx2,dy2,ds,l11,l12,l21,l22,sq1,sq2,eps=1.e-14;
    int j=0;
    jp = 1;
    for (int n=0; n<nFaces-nGfaces; n++){
      n1          = face(0,n);
      n2          = face(1,n);
      dx1         = x (0,j ,n2)-x (0,j,n1);
      dy1         = x (1,j ,n2)-x (1,j,n1);
      dx2         = xc(0,jp,n )-xc(0,j,n );
      dy2         = xc(1,jp,n )-xc(1,j,n );
      ds          = dx1*dy2-dx2*dy1;
      if (fabs(ds) < eps) ds = 0.; // on sharp corners qx = 0.
      else ds = 1./ds;
      l11         = ds*dy2;
      l12         =-ds*dy1;
      l21         =-ds*dx2;
      l22         = ds*dx1;
    for (int kk=0; kk<nqGradQ; kk++){
      k           = iqgrad(kk);
      sq1         = qp(k,j ,n2)-qp(k,j,n1);
      sq2         = q (k,jp,n )-q (k,j,n );
      qx(k,0,j,n) = l11*sq1+l12*sq2;
      qx(k,1,j,n) = l21*sq1+l22*sq2;
    }}

    j  = nPstr+1;
    jm = nPstr;
    for (int n=0; n<nFaces-nGfaces; n++){
      n1          = face(0,n);
      n2          = face(1,n);
      dx1         = x (0,j ,n2)-x (0,j,n1);
      dy1         = x (1,j ,n2)-x (1,j,n1);
      dx2         = xc(0,jm,n )-xc(0,j,n );
      dy2         = xc(1,jm,n )-xc(1,j,n );
      ds          = dx1*dy2-dx2*dy1;
      if (fabs(ds) < eps) ds = 0.;
      else ds = 1./ds;
      l11         = ds*dy2;
      l12         =-ds*dy1;
      l21         =-ds*dx2;
      l22         = ds*dx1;
    for (int kk=0; kk<nqGradQ; kk++){
      k           = iqgrad(kk);
      sq1         = qp(k,j ,n2)-qp(k,j,n1);
      sq2         = q (k,jm,n )-q (k,j,n );
      qx(k,0,j,n) = l11*sq1+l12*sq2;
      qx(k,1,j,n) = l21*sq1+l22*sq2;
    }}

    for (int n=nEdges-nBedges; n<nEdges; n++){
      c1           = edge(0,n);
      c2           = edge(1,n);
      n1           = edgn(  n);
    for (int j=1; j<nPstr+1; j++){
      jm           = j-1;
      dx1          = x (0,j,n1)-x (0,jm,n1);
      dy1          = x (1,j,n1)-x (1,jm,n1);
      dx2          = xc(0,j,c1)-xc(0,j ,c2);
      dy2          = xc(1,j,c1)-xc(1,j ,c2);
      ds  = dx1*dy2-dx2*dy1;
      if (fabs(ds) < eps) ds = 0.;
      else ds  = 1./ds;
      l11          = ds*dy2;
      l12          =-ds*dy1;
      l21          =-ds*dx2;
      l22          = ds*dx1;
    for (int kk=0; kk<nqGradQ; kk++){
      k            = iqgrad(kk);
      sq1          = qp(k,j,n1)-qp(k,jm,n1);
      sq2          = q (k,j,c1)-q (k,j ,c2);
      qx(k,0,j,c2) = l11*sq1+l12*sq2;
      qx(k,1,j,c2) = l21*sq1+l22*sq2;
    }}}
    gradQFlag = 1;
  }
}
