#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::faceAreaVis()
{
  // form face areas for viscous terms
  int n1,n2,n3,n4,j1,j2,j3,j4;
  double xa,ya,xb,yb,x1,y1,x2,y2,x3,y3,x4,y4;
  for (int n=0; n<nSurfElem; n++)
    for (int j=0; j<nStrandElem; j++)
      for (int i=0; i<nee; i++){
	n1           = edgeE(i,0,0);
	j1           = edgeE(i,0,1);
	n2           = edgeE(i,1,0);
	j2           = edgeE(i,1,1);
	n3           = edgeE(i,2,0);
	j3           = edgeE(i,2,1);
	n4           = edgeE(i,3,0);
	j4           = edgeE(i,3,1);

	n1           = surfElem(n,n1);
	j1           = j*meshOrder+j1;
	n2           = surfElem(n,n2);
	j2           = j*meshOrder+j2;
	x1           = x(strandMap(n1,j1),0);
	y1           = x(strandMap(n1,j1),1);
	x2           = x(strandMap(n2,j2),0);
	y2           = x(strandMap(n2,j2),1);

	if (n3 != -1){
	  n3         = surfElem(n,n3);
	  j3         = j*meshOrder+j3;
	  xa         = x(strandMap(n3,j3),0);
	  ya         = x(strandMap(n3,j3),1);
	}
	else{
	  xa         = .5*(x1+x2);
	  ya         = .5*(y1+y2);
	}

	if (n4 != -1){
	  n4         = surfElem(n,n4);
	  j4         = j*meshOrder+j4;
	  xb         = x(strandMap(n4,j4),0);
	  yb         = x(strandMap(n4,j4),1);
	}
	else{
	  xb         = .5*(x1+x2);
	  yb         = .5*(y1+y2);
	}

	areaE(n,j,i,0) =(ya-yb)/3.;
	areaE(n,j,i,1) =(xb-xa)/3.;
	//cout << n << " " << j << "   " << strandMap(n1,j1) << " " << strandMap(n2,j2) << " " << i << " " << areaE(n,j,i,0) << " " << areaE(n,j,i,1) << endl;
	}
}
