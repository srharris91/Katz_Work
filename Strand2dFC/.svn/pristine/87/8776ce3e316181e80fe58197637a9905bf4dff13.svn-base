#include "Strand2dFCBlockSolver.h"
#include "lagrangePoly.h"


void Strand2dFCBlockSolver::prepare()
{
  rms.allocate(nq);
  q.allocate(nSurfNode,nStrandNode,nq);
  qa.allocate(nSurfNode,nStrandNode,nqa);
  q0.allocate(nSurfNode,nStrandNode,nq);
  if (level > 0) fwc.allocate(nSurfNode,nStrandNode,nq);
  if (level > 0) fc.allocate(nSurfNode,nStrandNode,nq);
  qn.allocate(nSurfNode,nStrandNode,nq);
  if (nSteps > 0) qt.allocate(nSurfElem,meshOrder+1,nStrandNode,nq,timeAcc);
  dt.allocate(nSurfNode,nStrandNode);
  r.allocate(nSurfNode,nStrandNode,nq);
  r.set(0.);
  d.allocate(nSurfNode,nStrandNode,nq);
  d.set(0.);
  dn.allocate(nSurfNode,nStrandNode,nq);
  dn.set(0.);
  if (sourceMMS == 1) s.allocate(nSurfNode,nStrandNode,nq);
  qx.allocate(nSurfNode,nStrandNode,nq,2);
  lim.allocate(nSurfNode,nStrandNode,nq);
  surfData.allocate(nSurfNode,2,nq);
  if (viscous) surfDataVis.allocate(nSurfNode,2,nq,2);
  bndData.allocate(nBndNode,nStrandNode,nq);
  if (viscous) bndDataVis.allocate(nBndNode,nStrandNode,nq,2);
  Ads.allocate(nSurfElem,meshOrder+1,nq,nq);
  Adl.allocate(nSurfNode,nq,nq);
  if (viscous) qas.allocate(nSurfElem,meshOrder+1,nStrandNode,nqaGradQa);


  // read q data from restart
  if (restartStep > 0){
    cout << "\n*** restart capability not yet implemented in prepare.C ***"
         << endl;
    exit(0);
  }


  // initialize q and qa at interior and flux point locations
  double x0,y0,nx,ny;
  Array3D<double> x(nSurfNode,nStrandNode,2);
  Array2D<double> dw(nSurfNode,nStrandNode);
  for (int n=0; n<nSurfNode; n++){
    x0 = surfX(n,0);
    y0 = surfX(n,1);
    nx = pointingVec(n,0);
    ny = pointingVec(n,1);
    for (int j=0; j<nStrandNode; j++){
      x(n,j,0) = x0+nx*strandX(j);
      x(n,j,1) = y0+ny*strandX(j);
      dw(n,j)  = strandX(j); //strand length
    }}
  //if (sourceMMS == 1) dw.set(1.e14);
  sys->initWallDist(nSurfNode*nStrandNode,
                    &dw(0,0),
                    &qa(0,0,0));
  sys->initFlow(nSurfNode*nStrandNode,
                &x(0,0,0),
                &q(0,0,0));
  sys->stepQAdd(nSurfNode*nStrandNode,
                &q(0,0,0),
                &qa(0,0,0));


  // initialize MMS source terms at interior nodes
  if (sourceMMS == 1){
    s.set(0.);
    sys->initSource(nSurfNode*nStrandNode,
                    &x(0,0,0),
                    &s(0,0,0));
  }


  // initialize boundary data to the exact MMS solution
  if (sourceMMS == 1){
    for (int n=0; n<nSurfNode; n++){
      for (int k=0; k<nq; k++) surfData(n,0,k) = q(n,0            ,k);
      for (int k=0; k<nq; k++) surfData(n,1,k) = q(n,nStrandNode-1,k);
    }
    int m;
    for (int n=0; n<nBndNode; n++){
      m = bndNode(n);
      for (int j=0; j<nStrandNode; j++)
	for (int k=0; k<nq; k++) bndData(n,j,k) = q(m,j,k);
    }
    if (viscous){
      for (int n=0; n<nSurfNode; n++){
	sys->initVisFlux(1,&x(n,0            ,0),&surfDataVis(n,0,0,0));
	sys->initVisFlux(1,&x(n,nStrandNode-1,0),&surfDataVis(n,1,0,0));
      }
      for (int n=0; n<nBndNode; n++){
	m = bndNode(n);
	for (int j=0; j<nStrandNode; j++)
	  sys->initVisFlux(1,&x(m,j,0),&bndDataVis(n,j,0,0));
      }}
  }

  else{ // set boundary data using boundary condition data in input file
    for (int n=0; n<nSurfNode; n++){
      sys->initPenaltyData(1,&surfNodeTag(n,0),
			   &x(n,0            ,0),&surfData(n,0,0));
      sys->initPenaltyData(1,&surfNodeTag(n,1),
			   &x(n,nStrandNode-1,0),&surfData(n,1,0));
    }
    int m;
    for (int n=0; n<nBndNode; n++){
      m = bndNode(n);
      for (int j=0; j<nStrandNode; j++)
	sys->initPenaltyData(1,&bndNodeTag(n),&x(m,j,0),&bndData(n,j,0));
    }
    if (viscous){
      surfDataVis.set(0.); //assume far field conditions have 0 viscous flux
      bndDataVis.set(0.);
    }
  }

  x.deallocate();
 dw.deallocate();
}
