#include "StrandBlockSolver.h"


void StrandBlockSolver::finalize(const int& mglevel)
{
  if (sys){
    if (mglevel == 0){
      sys->finalize();
      delete sys;
    }
    sys = NULL;
  }
  iqgrad.deallocate();
  iqagrad.deallocate();
  rmsNorm.deallocate();
  rms.deallocate();
  face.deallocate();
  fTag.deallocate();
  bTag.deallocate();
  fClip.deallocate();
  sFlag.deallocate();
  edge.deallocate();
  edgp.deallocate();
  edgn.deallocate();
  ncsc.deallocate();
  csc.deallocate();
  ncsp.deallocate();
  if (csp){
    for (int n=0; n<nNodes; n++) if (csp[n]) delete [] csp[n];
    delete [] csp;
    csp = NULL;
  }
  gsMap.deallocate();
  f2cc.deallocate();
  f2ce.deallocate();
  f2cs.deallocate();
  x.deallocate();
  x0.deallocate();
  x1.deallocate();
  x2.deallocate();
  xc.deallocate();
  v.deallocate();
  v1.deallocate();
  v2.deallocate();
  facu.deallocate();
  facs.deallocate();
  if (lsp){
    for (int n=0; n<2*(nPstr+1)*nNodes; n++) if (lsp[n]) delete [] lsp[n];
    delete [] lsp;
    lsp = NULL;
  }
  xStr.deallocate();
  xvu.deallocate();
  xvs.deallocate();
  dvu.deallocate();
  dvs.deallocate();
  nvu.deallocate();
  nvs.deallocate();
  q.deallocate();
  qa.deallocate();
  q0.deallocate();
  q1.deallocate();
  q2.deallocate();
  qp.deallocate();
  qap.deallocate();
  qx.deallocate();
  qax.deallocate();
  r.deallocate();
  dq.deallocate();
  s.deallocate();
  dt.deallocate();
  radi.deallocate();
  radv.deallocate();
  rads.deallocate();
  dd.deallocate();
  dm.deallocate();
  dp.deallocate();
  bu.deallocate();
  fwc.deallocate();
  limu.deallocate();
  lims.deallocate();
  dlim.deallocate();

  dataInit();
}
