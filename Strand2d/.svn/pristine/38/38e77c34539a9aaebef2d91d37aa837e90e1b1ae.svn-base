#include "StrandBlockSolver.h"


void StrandBlockSolver::gradQa(const int& mglevel)
{
  if (mglevel == 0 && gradient != 0 && gradQaFlag == 0){

    // compute nodal values
    nodalQa(mglevel);


    // initialize gradients
    for (int n=0; n<nFaces+nBedges; n++)
      for (int j=0; j<nPstr+2; j++)
	for (int k=0; k<ndim; k++)
	  for (int m=0; m<nqa; m++) qax(m,k,j,n) = 0.;


    // loop through unstructured edges
    int c1,c2,n1,n2,jm,jp,k;
    double Ax,Ay,sqa;
    for (int n=0; n<nEdges; n++){
      c1             = edge(0,n);
      c2             = edge(1,n);
      n1             = edgn(n);
    for (int j=1; j<nPstr+1; j++){
      jm             = j-1;
      Ax             = facs(0,j,n);
      Ay             = facs(1,j,n);
    for (int kk=0; kk<nqaGradQa; kk++){
      k              = iqagrad(kk);
      sqa            = qap(k,jm,n1)+qap(k,j,n1);
      qax(k,0,j,c1) += Ax*sqa;
      qax(k,1,j,c1) += Ay*sqa;
      qax(k,0,j,c2) -= Ax*sqa;
      qax(k,1,j,c2) -= Ay*sqa;
    }}}


    // loop through structured edges
    for (int n=0; n<nFaces-nGfaces; n++){
      n1             = face(0,n);
      n2             = face(1,n);
    for (int j=0; j<nPstr+1; j++){
      jp             = j+1;
      Ax             = facu(0,j,n);
      Ay             = facu(1,j,n);
    for (int kk=0; kk<nqaGradQa; kk++){
      k              = iqagrad(kk);
      sqa            = qap(k,j,n1)+qap(k,j,n2);
      qax(k,0,j ,n) += Ax*sqa;
      qax(k,1,j ,n) += Ay*sqa;
      qax(k,0,jp,n) -= Ax*sqa;
      qax(k,1,jp,n) -= Ay*sqa;
    }}}


    // divide by twice the volume for the interior cells
    for (int n=0; n<nFaces-nGfaces; n++)
      for (int j=1; j<nPstr+1; j++)
	for (int m=0; m<ndim; m++)
	  for (int kk=0; kk<nqaGradQa; kk++){
	    k             = iqagrad(kk);
	    qax(k,m,j,n) /= (2.*v(j,n));
	  }


    // surface, end, and boundary gradients
    double dx1,dy1,dx2,dy2,ds,l11,l12,l21,l22,sqa1,sqa2,eps=1.e-14;
    int j=0;
    jp = 1;
    for (int n=0; n<nFaces-nGfaces; n++){
      n1           = face(0,n);
      n2           = face(1,n);
      dx1          = x (0,j ,n2)-x (0,j,n1);
      dy1          = x (1,j ,n2)-x (1,j,n1);
      dx2          = xc(0,jp,n )-xc(0,j,n );
      dy2          = xc(1,jp,n )-xc(1,j,n );
      ds           = dx1*dy2-dx2*dy1;
      if (fabs(ds) < eps) ds = 0.; // on sharp corners qax = 0.
      else ds = 1./ds;
      l11          = ds*dy2;
      l12          =-ds*dy1;
      l21          =-ds*dx2;
      l22          = ds*dx1;
    for (int kk=0; kk<nqaGradQa; kk++){
      k            = iqagrad(kk);
      sqa1         = qap(k,j ,n2)-qap(k,j,n1);
      sqa2         = qa (k,jp,n )-qa (k,j,n );
      qax(k,0,j,n) = l11*sqa1+l12*sqa2;
      qax(k,1,j,n) = l21*sqa1+l22*sqa2;
    }}

    j  = nPstr+1;
    jm = nPstr;
    for (int n=0; n<nFaces-nGfaces; n++){
      n1           = face(0,n);
      n2           = face(1,n);
      dx1          = x (0,j ,n2)-x (0,j,n1);
      dy1          = x (1,j ,n2)-x (1,j,n1);
      dx2          = xc(0,jm,n )-xc(0,j,n );
      dy2          = xc(1,jm,n )-xc(1,j,n );
      ds           = dx1*dy2-dx2*dy1;
      if (fabs(ds) < eps) ds = 0.;
      else ds = 1./ds;
      l11          = ds*dy2;
      l12          =-ds*dy1;
      l21          =-ds*dx2;
      l22          = ds*dx1;
    for (int kk=0; kk<nqaGradQa; kk++){
      k            = iqagrad(kk);
      sqa1         = qap(k,j ,n2)-qap(k,j,n1);
      sqa2         = qa (k,jm,n )-qa (k,j,n );
      qax(k,0,j,n) = l11*sqa1+l12*sqa2;
      qax(k,1,j,n) = l21*sqa1+l22*sqa2;
    }}

    for (int n=nEdges-nBedges; n<nEdges; n++){
      c1            = edge(0,n);
      c2            = edge(1,n);
      n1            = edgn(  n);
    for (int j=1; j<nPstr+1; j++){
      jm            = j-1;
      dx1           = x (0,j,n1)-x (0,jm,n1);
      dy1           = x (1,j,n1)-x (1,jm,n1);
      dx2           = xc(0,j,c1)-xc(0,j ,c2);
      dy2           = xc(1,j,c1)-xc(1,j ,c2);
      ds  = dx1*dy2-dx2*dy1;
      if (fabs(ds) < eps) ds = 0.;
      else ds  = 1./ds;
      l11           = ds*dy2;
      l12           =-ds*dy1;
      l21           =-ds*dx2;
      l22           = ds*dx1;
    for (int kk=0; kk<nqaGradQa; kk++){
      k             = iqagrad(kk);
      sqa1          = qap(k,j,n1)-qap(k,jm,n1);
      sqa2          = qa (k,j,c1)-qa (k,j ,c2);
      qax(k,0,j,c2) = l11*sqa1+l12*sqa2;
      qax(k,1,j,c2) = l21*sqa1+l22*sqa2;
    }}}
    gradQaFlag = 1;
  }
}
