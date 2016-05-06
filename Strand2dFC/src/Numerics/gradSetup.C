#include "Strand2dFCBlockSolver.h"
#include "lagrangePoly.h"


void Strand2dFCBlockSolver::gradSetup()
{
  gxc.allocate(npsp1,nStrandNode,2);

  // use linear gradients if second-order in the surface direction
  if (surfOrder == 2){
    Array3D<double> gxL(npsp1,nStrandNode,2);
    gradSetupFull(gxL);
    for (int n=0; n<npsp1; n++)
      for (int j=0; j<nStrandNode; j++)
	for (int k=0; k<2; k++) gxc(n,j,k) = gxL(n,j,k);
    gxL.deallocate();
  }


  // use high order mapping if high-order in the surface direction
  else{
    Array3D<double> gxQ(npsp1,nStrandNode,2);
    Array3D<double> gxC(npsp1,nStrandNode,2);

    // set up averaged quadratic gradient coefficients
    if (gradMethod == 0 || gradMethod == 2) gradSetupSub(gxQ);


    // set up full order gradient coefficients
    if (gradMethod == 1 || gradMethod == 2) gradSetupFull(gxC);


    // copy the proper coefficients based on desired gradient method
    int m;
    if      (gradMethod == 0)
      for (int n=0; n<npsp1; n++)
	for (int j=0; j<nStrandNode; j++)
	  for (int k=0; k<2; k++) gxc(n,j,k) = gxQ(n,j,k);
    else if (gradMethod == 1)
      for (int n=0; n<npsp1; n++)
	for (int j=0; j<nStrandNode; j++)
	  for (int k=0; k<2; k++) gxc(n,j,k) = gxC(n,j,k);
    else if (gradMethod == 2){
      for (int n=0; n<npsp1; n++)
	for (int j=0; j<nStrandNode; j++)
	  for (int k=0; k<2; k++) gxc(n,j,k) = gxC(n,j,k);
      for (int n=0; n<nBndNode; n++){
	m = bndNode(n);
	for(int i=psp2(m); i<psp2(m+1); i++)
	  for (int j=0; j<nStrandNode; j++)
	    for (int k=0; k<2; k++) gxc(i,j,k) = gxQ(i,j,k);
      }}
    else{
      cout << "\ngradMethod not recognized in gradSetup.C" << endl;
      exit(0);
    }
    gxQ.deallocate();
    gxC.deallocate();
  }

  /*
  for (int n=0; n<nSurfNode; n++)
    for (int i=psp2(n); i<psp2(n+1); i++)
      cout << n << " " << psp1(i) << " " << gxc(i,0,0) << " " << gxc(i,0,1) << endl;
  exit(0);
  */
}
