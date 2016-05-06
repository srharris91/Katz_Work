#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::finalize()
{
  orders.deallocate();

  if (sys){
    sys->finalize();
    delete sys;
    sys = NULL;
  }

  if (triG){
    for (int n=0; n<nTriG; n++){
      if (triG[n]) delete [] triG[n];
      triG[n] = NULL;
    }
    delete [] triG;
    triG = NULL;
  }
  xG.deallocate();
  edgeBdG.deallocate();

  if (nfn) delete [] nfn;
  nfn = NULL;
  if (nfn1){
    for (int n=0; n<nNodeF; n++){
      if (nfn1[n]) delete [] nfn1[n];
      nfn1[n] = NULL;
    }
    delete [] nfn1;
    nfn1 = NULL;
  }
  if (nfn2){
    for (int n=0; n<nNodeF; n++){
      if (nfn2[n]) delete [] nfn2[n];
      nfn2[n] = NULL;
    }
    delete [] nfn2;
    nfn2 = NULL;
  }
  lqCF.deallocate();
  lqFC.deallocate();
  nfe.deallocate();
  nce.deallocate();

  iqgrad.deallocate();
  iqagrad.deallocate();
  elem.deallocate();
  tri.deallocate();
  edge.deallocate();
  edgeE.deallocate();
  edgeBd.deallocate();
  nodeBd.deallocate();
  edgeQ.deallocate();
  psp1.deallocate();
  psp2.deallocate();
  x.deallocate();
  v.deallocate();
  area.deallocate();
  areaBd.deallocate();
  areaE.deallocate();
  rka.deallocate();
  rkb.deallocate();
  bdf.deallocate();
  dlim.deallocate();
  rmsNorm.deallocate();
  rms.deallocate();
  ln.deallocate();
  gxS.deallocate();
  nxQ.deallocate();
  dxg.deallocate();
  gx.deallocate();
  gxQ.deallocate();
  gxC.deallocate();
  gxx.deallocate();
  outputVarLength.deallocate();
  outputVars.deallocate();
  q.deallocate();
  qa.deallocate();
  q0.deallocate();
  fwc.deallocate();
  qn.deallocate();
  qt.deallocate();
  radi.deallocate();
  radv.deallocate();
  dt.deallocate();
  r.deallocate();
  d.deallocate();
  dn.deallocate();
  s.deallocate();
  qx.deallocate();
  lim.deallocate();

  dataInit();
}
