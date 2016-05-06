#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::pseudoTime(const int& j)
{
  for (int n=0; n<nSurfNode; n++) dt(n,j) = 0.;


  // inviscid spectral radius
  if (inviscid == 1){
    radi.allocate(nSurfNode);
    specRadi(j);
    for (int n=0; n<nSurfNode; n++) dt(n,j) = radi(n)/cfl;
    radi.deallocate();
  }


  // viscous spectral radius
  if (viscous == 1){
    radv.allocate(nSurfNode);
    specRadv(j);
    for (int n=0; n<nSurfNode; n++) dt(n,j) = max(radv(n)/vnn,dt(n,j));
    radv.deallocate();
  }


  // pseudo-time step estimate
  for (int n=0; n<nSurfNode; n++) dt(n,j) = v(n,j)/dt(n,j);
  //double a=1.;
  //for(int n=0; n<nSurfNode; n++) if (dt(n,j) < a) a = dt(n,j);
  //cout << j << " " << a << endl;
}
