#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::pseudoTime()
{
  dt.set(0.);

  // inviscid spectral radius
  if (inviscid == 1){
    radi.allocate(nNode);
    specRadi();
    for(int n=0; n<nNode; n++) dt(n) = radi(n)/cfl;
    radi.deallocate();
  }


  // viscous spectral radius
  if (viscous == 1){
    radv.allocate(nNode);
    specRadv();
    for(int n=0; n<nNode; n++) dt(n) = max(radv(n)/vnn,dt(n));
    radv.deallocate();
  }


  // pseudo-time step estimate
  for(int n=0; n<nNode; n++) dt(n) = v(n)/dt(n);


  /*
  // set the time step to the minimum time step over the mesh
  double dtMin=dt(0);
  for (int n=0; n<nNode; n++) dtMin = min(dtMin,dt(n));
  dt.set(dtMin);
  */
}
