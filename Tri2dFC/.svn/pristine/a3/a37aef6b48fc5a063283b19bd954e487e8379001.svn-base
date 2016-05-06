#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::faceAreaVis()
{
  // form face areas for viscous terms
  int k,orderE=3-level; //order of local elements
  Array2D<double> xE(nne,2);
  double xa,ya,xb,yb;
  if      (orderE == 1)
    for (int n=0; n<nElem; n++){
      for (int j=0; j<nne; j++){
	xE(j,0) = x(elem(n,j),0);
	xE(j,1) = x(elem(n,j),1);
      }
      k              = 0;
      xa             = .5*(xE(0,0)+xE(1,0));
      ya             = .5*(xE(0,1)+xE(1,1));
      xb             = xE(2,0);
      yb             = xE(2,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = xE(1,0);
      ya             = xE(1,1);
      xb             = .5*(xE(0,0)+xE(2,0));
      yb             = .5*(xE(0,1)+xE(2,1));
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = .5*(xE(1,0)+xE(2,0));
      ya             = .5*(xE(1,1)+xE(2,1));
      xb             = xE(0,0);
      yb             = xE(0,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
    }
  
  else if (orderE == 2)
    for (int n=0; n<nElem; n++){
      for (int j=0; j<nne; j++){
	xE(j,0) = x(elem(n,j),0);
	xE(j,1) = x(elem(n,j),1);
      }
      k              = 0;
      xa             = .5*(xE(0,0)+xE(3,0));
      ya             = .5*(xE(0,1)+xE(3,1));
      xb             = xE(5,0);
      yb             = xE(5,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = .5*(xE(3,0)+xE(1,0));
      ya             = .5*(xE(3,1)+xE(1,1));
      xb             = xE(4,0);
      yb             = xE(4,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = xE(3,0);
      ya             = xE(3,1);
      xb             = .5*(xE(0,0)+xE(5,0));
      yb             = .5*(xE(0,1)+xE(5,1));
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = xE(4,0);
      ya             = xE(4,1);
      xb             = xE(0,0);
      yb             = xE(0,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = xE(1,0);
      ya             = xE(1,1);
      xb             = xE(5,0);
      yb             = xE(5,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = .5*(xE(1,0)+xE(4,0));
      ya             = .5*(xE(1,1)+xE(4,1));
      xb             = xE(3,0);
      yb             = xE(3,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = xE(3,0);
      ya             = xE(3,1);
      xb             = xE(2,0);
      yb             = xE(2,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = xE(4,0);
      ya             = xE(4,1);
      xb             = .5*(xE(2,0)+xE(5,0));
      yb             = .5*(xE(2,1)+xE(5,1));
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = .5*(xE(2,0)+xE(4,0));
      ya             = .5*(xE(2,1)+xE(4,1));
      xb             = xE(5,0);
      yb             = xE(5,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
    }
  
  else if (orderE == 3)
    for (int n=0; n<nElem; n++){
      for (int j=0; j<nne; j++){
	xE(j,0) = x(elem(n,j),0);
	xE(j,1) = x(elem(n,j),1);
      }
      k              = 0;
      xa             = .5*(xE(0,0)+xE(3,0));
      ya             = .5*(xE(0,1)+xE(3,1));
      xb             = xE(8,0);
      yb             = xE(8,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = .5*(xE(3,0)+xE(4,0));
      ya             = .5*(xE(3,1)+xE(4,1));
      xb             = xE(9,0);
      yb             = xE(9,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = .5*(xE(4,0)+xE(1,0));
      ya             = .5*(xE(4,1)+xE(1,1));
      xb             = xE(5,0);
      yb             = xE(5,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = xE(3,0);
      ya             = xE(3,1);
      xb             = .5*(xE(0,0)+xE(8,0));
      yb             = .5*(xE(0,1)+xE(8,1));
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = xE(9,0);
      ya             = xE(9,1);
      xb             = xE(0,0);
      yb             = xE(0,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = xE(4,0);
      ya             = xE(4,1);
      xb             = xE(8,0);
      yb             = xE(8,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = xE(5,0);
      ya             = xE(5,1);
      xb             = xE(3,0);
      yb             = xE(3,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = xE(1,0);
      ya             = xE(1,1);
      xb             = xE(9,0);
      yb             = xE(9,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = .5*(xE(1,0)+xE(5,0));
      ya             = .5*(xE(1,1)+xE(5,1));
      xb             = xE(4,0);
      yb             = xE(4,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = xE(3,0);
      ya             = xE(3,1);
      xb             = xE(7,0);
      yb             = xE(7,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = xE(4,0);
      ya             = xE(4,1);
      xb             = xE(6,0);
      yb             = xE(6,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = xE(9,0);
      ya             = xE(9,1);
      xb             = .5*(xE(7,0)+xE(8,0));
      yb             = .5*(xE(7,1)+xE(8,1));
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = xE(6,0);
      ya             = xE(6,1);
      xb             = xE(8,0);
      yb             = xE(8,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = xE(5,0);
      ya             = xE(5,1);
      xb             = xE(7,0);
      yb             = xE(7,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = .5*(xE(5,0)+xE(6,0));
      ya             = .5*(xE(5,1)+xE(6,1));
      xb             = xE(9,0);
      yb             = xE(9,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = xE(9,0);
      ya             = xE(9,1);
      xb             = xE(2,0);
      yb             = xE(2,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = xE(6,0);
      ya             = xE(6,1);
      xb             = .5*(xE(2,0)+xE(7,0));
      yb             = .5*(xE(2,1)+xE(7,1));
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
      xa             = .5*(xE(6,0)+xE(2,0));
      ya             = .5*(xE(6,1)+xE(2,1));
      xb             = xE(7,0);
      yb             = xE(7,1);
      areaE(n,k  ,0) = yb-ya;
      areaE(n,k++,1) = xa-xb;
    }
  
  for (int n=0; n<nElem; n++)
    for (int j=0; j<nee; j++){
      areaE(n,j,0) /= 3.;
      areaE(n,j,1) /= 3.;
    }
  xE.deallocate();
}
