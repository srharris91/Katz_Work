#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::gradient(Array3D<double>& p,
				     Array4D<double>& px,
				     const int& j)
{
  for (int n=0; n<nSurfNode; n++)
    for (int k=0; k<nq; k++){
      px(n,j,k,0) = 0.;
      px(n,j,k,1) = 0.;
    }

  if (surfOrder > 1){
    int ni;
    double cx,cy;
    for (int n=0; n<nSurfNode; n++)
      for (int i=psp2(n); i<psp2(n+1); i++){
	ni = psp1(i);
	cx = gxc(i,j,0);
	cy = gxc(i,j,1);
	for (int k=0; k<nq; k++){	
	  px(n,j,k,0) += cx*p(ni,j,k);
	  px(n,j,k,1) += cy*p(ni,j,k);
	}}}
}
