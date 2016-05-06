#include "Tri2dFCBlockSolver.h"
#include "lagrangePoly.h"


void Tri2dFCBlockSolver::gradSetup()
{
  // set up averaged quadratic gradient coefficients
  if (gradMethod == 0 || gradMethod == 2){
    gxQ.allocate(npsp1,2);
    gradSetupQuadratic();
  }


  // set up fully cubic gradient coefficients
  gxC.allocate(npsp1,2);
  gradSetupCubic();


  // copy the proper coefficients based on desired gradient method
  if      (gradMethod == 0)
    for (int n=0; n<npsp1; n++)
      for (int k=0; k<2; k++) gx(n,k) = gxQ(n,k);
  else if (gradMethod == 1)
    for (int n=0; n<npsp1; n++)
      for (int k=0; k<2; k++) gx(n,k) = gxC(n,k);
  else if (gradMethod == 2){
    for (int n=0; n<npsp1; n++)
      for (int k=0; k<2; k++) gx(n,k) = gxC(n,k);
    for (int n=nNode-nNodeBd; n<nNode; n++)
      for(int i=psp2(n); i<psp2(n+1); i++)
	for (int k=0; k<2; k++) gx(i,k) = gxQ(i,k);
  }
  else{
    cout << "\ngradMethod not recognized in gradSetup.C" << endl;
    exit(0);
  }


  // deallocate work arrays
  gxQ.deallocate();
  gxC.deallocate();
}
