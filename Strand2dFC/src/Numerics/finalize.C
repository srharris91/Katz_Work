#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::finalize()
{
  if (sys){
    sys->finalize();
    delete sys;
    sys = NULL;
  }

  if (icn1){
    for (int j=0; j<nStrandNode; j++){
      if (icn1[j]) delete [] icn1[j];
      icn1[j] = NULL;
    }
    delete [] icn1;
    icn1 = NULL;
  }
  if (icn2){
    for (int j=0; j<nStrandNode; j++){
      if (icn2[j]) delete [] icn2[j];
      icn2[j] = NULL;
    }
    delete [] icn2;
    icn2 = NULL;
  }
  if (dcn1){
    for (int j=0; j<nStrandNode+1; j++){
      if (dcn1[j]) delete [] dcn1[j];
      dcn1[j] = NULL;
    }
    delete [] dcn1;
    dcn1 = NULL;
  }
  if (dcn2){
    for (int j=0; j<nStrandNode+1; j++){
      if (dcn2[j]) delete [] dcn2[j];
      dcn2[j] = NULL;
    }
    delete [] dcn2;
    dcn2 = NULL;
  }

  if (vcn1){
    for (int j=0; j<nStrandNode; j++){
      if (vcn1[j]){
	int nj=vcn2[j][1];
	for (int m=0; m<nj; m++){
	  if (vcn1[j][m]) delete [] vcn1[j][m];
	  vcn1[j][m] = NULL;
	}
	delete [] vcn1[j];
	vcn1[j] = NULL;
      }}
    delete [] vcn1;
    vcn1 = NULL;
  }
  if (vcn3){
    for (int j=0; j<nStrandNode; j++){
      if (vcn3[j]){
	int nj=vcn2[j][1];
	for (int m=0; m<nj; m++){
	  if (vcn3[j][m]) delete [] vcn3[j][m];
	  vcn3[j][m] = NULL;
	}
	delete [] vcn3[j];
	vcn3[j] = NULL;
      }}
    delete [] vcn3;
    vcn3 = NULL;
  }
  if (vcn2){
    for (int j=0; j<nStrandNode; j++){
      if (vcn2[j]) delete [] vcn2[j];
      vcn2[j] = NULL;
    }
    delete [] vcn2;
    vcn2 = NULL;
  }
  if (vcn4) delete [] vcn4;
  vcn4 = NULL;

  surfElem.deallocate();
  surfElem0.deallocate();
  surfElemTag.deallocate();
  surfNodeTag.deallocate();
  bndNode.deallocate();
  bndElem.deallocate();
  bndNodeTag.deallocate();
  bndNodeNormal.deallocate();
  surfX.deallocate();
  strandX.deallocate();
  pointingVec.deallocate();
  clip.deallocate();
  iqgrad.deallocate();
  iqagrad.deallocate();
  rka.deallocate();
  rkb.deallocate();
  bdf.deallocate();
  dlim.deallocate();
  rmsNorm.deallocate();
  rms.deallocate();
  outputVarLength.deallocate();
  outputVars.deallocate();
  elemEdge.deallocate();
  surfEdge.deallocate();
  bndSign.deallocate();
  psp1.deallocate();
  psp2.deallocate();
  psp1S.deallocate();
  psp2S.deallocate();
  wsp1S.deallocate();
  psp1MGS.deallocate();
  psp2MGS.deallocate();
  wsp1MGS.deallocate();
  psp1MGR.deallocate();
  psp2MGR.deallocate();
  wsp1MGR.deallocate();
  psp1MGC.deallocate();
  psp2MGC.deallocate();
  wsp1MGC.deallocate();
  ls.deallocate();
  xn.deallocate();
  yn.deallocate();
  xs.deallocate();
  ys.deallocate();
  xsA.deallocate();
  ysA.deallocate();
  xns.deallocate();
  yns.deallocate();
  jac.deallocate();
  v.deallocate();
  sn.deallocate();
  wQ.deallocate();
  lQ.deallocate();
  lsQ.deallocate();
  xsQ.deallocate();
  ysQ.deallocate();
  xnQ.deallocate();
  ynQ.deallocate();
  jacQ.deallocate();
  gxc.deallocate();
  q.deallocate();
  qa.deallocate();
  q0.deallocate();
  fc.deallocate();
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
  surfData.deallocate();
  surfDataVis.deallocate();
  bndData.deallocate();
  bndDataVis.deallocate();
  Ads.deallocate();
  Adl.deallocate();
  qas.deallocate();

  dataInit();
}
