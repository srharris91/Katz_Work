#include "Strand2dFCBlockSolver.h"
#include "Strand2dFCSPLam.h"
#include "Strand2dFCSPTurbSA.h"
#include "Strand2dFCInputRead.h"


void Strand2dFCBlockSolver::initialize(const string& inputFile,
				       const int& level0,
				       const int& meshOrder00,
				       const int& nSurfElem00,
				       const int& nSurfNode0,
				       const int& nBndNode0,
				       const int& nStrandNode0,
				       const int& nCompBd0,
				       const int* surfElem00,
				       const int* bndNode0,
				       const double* surfX0,
				       const double* strandX0,
				       const int* surfElemTag0,
				       const int* bndNodeTag0,
				       const double* bndNodeNormal0,
				       const double* pointingVec0,
				       const int* clip0)
{
  // copy dimensions from mesh manager
  level       = level0;
  meshOrder0  = meshOrder00;
  nSurfElem0  = nSurfElem00;
  meshOrder   = meshOrder0; //assume this for now
  nSurfElem   = nSurfElem0; //assume this for now
  nSurfNode   = nSurfNode0;
  nStrandNode = nStrandNode0;
  nBndNode    = nBndNode0;
  nDofs       = nSurfNode*nStrandNode;


  // read input file
  int len=80;
  char systemType[len];
  int surfOrders[len];
  int strandOrders[len];
  strand2dfcinputread_(inputFile.size(),
		       inputFile.c_str(),
		       len,
		       iPrint,
		       iTest,
		       iDebug,
		       iConvFile,
		       iSolnFile,
		       iResdFile,
		       iErrFile,
		       iSurfFile,
		       standAlone,
		       restartStep,
		       nRestart,
		       nOutput,
		       nSteps,
		       nPseudoSteps,
		       nPseudoSteps0,
		       nLinearSteps,
		       nRKStages,
		       implicit,
		       &surfOrders[0],
		       &strandOrders[0],
		       gradMethod,
		       mgCycle,
		       limiter,
		       timeAcc,
		       dtUnsteady,
		       cfl,
		       vnn,
		       smooth,
		       convLimit,
		       relax,
		       systemType);

  surfOrder   = surfOrders[level];
  strandOrder = strandOrders[level];


  // allocate system layer instance and read system layer inputs
  if (!strcmp(systemType,"SPLam"))
    sys = new Strand2dFCSPLam;
  else if (!strcmp(systemType,"SPTurbSA"))
    sys = new Strand2dFCSPTurbSA;
  else{
    cout << "\nSystem layer not recognized in input.C: Terminating." << endl;
    exit(0);
  }


  sys->inputRead(inputFile);


  // setup the system layer
  int    tmp = 500;
  int    iqgradT [tmp];
  int    iqagradT[tmp];
  double rmsNormT[tmp];
  double dlimT   [tmp];
  int    outputVarLengthT[tmp];
  string outputVarsT[500];
  sys->prepSetup(iPrint,
                 iTest,
                 iDebug,
                 tmp,
		 iSolnFile,
		 iResdFile,
		 iErrFile,
                 nq,
                 nqa,
                 ndim,
                 inviscid,
                 viscous,
                 source,
                 sourceMMS,
                 dissipation,
                 nCompBd,
                 &iqgradT[0],
                 &iqagradT[0],
                 &dlimT[0],
                 &rmsNormT[0],
		 nOutputVars,
		 &outputVarLengthT[0],
		 &outputVarsT[0]);

  nOutputScalars = 0;
  nOutputVectors = 0;
  outputVarLength.allocate(nOutputVars);
  outputVars.allocate(nOutputVars);
  for (int n=0; n<nOutputVars; n++){
    outputVarLength(n) = outputVarLengthT[n];
    outputVars(n) = outputVarsT[n];
    if      (outputVarLength(n) == 1) nOutputScalars++;
    else if (outputVarLength(n) == 3) nOutputVectors++;
    else{
      cout << "\n*** Incorrect outputVarLength in input.C ***" << endl;
      exit(0);
    }}

  nqGradQ = 0;
  for (int n=0; n<nq; n++) if (iqgradT[n] == 1) nqGradQ++;
  iqgrad.allocate(nqGradQ);
  nqGradQ = 0;
  for (int n=0; n<nq; n++)
    if (iqgradT[n] == 1) iqgrad(nqGradQ++) = n;

  nqaGradQa = 0;
  for (int n=0; n<nqa; n++) if (iqagradT[n] == 1) nqaGradQa++;
  iqagrad.allocate(nqaGradQa);
  nqaGradQa = 0;
  for (int n=0; n<nqa; n++)
    if (iqagradT[n] == 1) iqagrad(nqaGradQa++) = n;

  rmsNorm.allocate(nq);
  for (int n=0; n<nq; n++) rmsNorm(n) = rmsNormT[n];

  dlim.allocate(nq);
  for (int n=0; n<nq; n++) dlim(n) = dlimT[n];

  rka.allocate(nRKStages);
  rkb.allocate(nRKStages);
  if (nRKStages == 5){
    rka(0) = .25;
    rka(1) = 1./6.;
    rka(2) = .375;
    rka(3) = .5;
    rka(4) = 1.;
    rkb(0) = 1.;
    rkb(1) = 0.;
    rkb(2) = .56;
    rkb(3) = 0.;
    rkb(4) = .44;
  }
  else{
    //cout << "\nCan only handle nRKStages = 5. Terminating" << endl;
    //exit(0);
    rka(0) = .25;
    rkb(0) = 1.;
  }

  bdf.allocate(timeAcc+1);
  if      (timeAcc == 1){
    bdf(0) = 1.;
    bdf(1) =-1.;
  }
  else if (timeAcc == 2){
    bdf(0) = 1.5;
    bdf(1) =-2.;
    bdf(2) = .5;
  }
  else if (timeAcc == 3){
    bdf(0) = 11./6.;
    bdf(1) =-3.;
    bdf(2) = 1.5;
    bdf(3) =-1./3.;
  }
  else{
    cout << "\nCan only handle timeAcc = 1, 2, or 3. Terminating" << endl;
    exit(0);
  }


  // copy needed data from the mesh manager
  if (nCompBd0+1 != nCompBd){ // I added strand ends
    cout << "\nNumber of boundary patches in the input file not consistent "
	 << "with the mesh file.\n" << endl;
    exit(0);
  }
  surfElem.allocate(nSurfElem,meshOrder+1);
  surfElemTag.allocate(nSurfElem);
  for (int n=0; n<nSurfElem; n++){
    for (int i=0; i<meshOrder+1; i++)
      surfElem(n,i) = surfElem00[n*(meshOrder+1)+i];
    surfElemTag(n) = surfElemTag0[n];
  }
  surfElem0.allocate(nSurfElem0,meshOrder0+1);
  for (int n=0; n<nSurfElem0; n++)
    for (int i=0; i<meshOrder0+1; i++) surfElem0(n,i) = surfElem(n,i);
  bndNode.allocate(nBndNode);
  bndNodeTag.allocate(nBndNode);
  bndNodeNormal.allocate(nBndNode,2);
  for (int n=0; n<nBndNode; n++){
    bndNode(n)         = bndNode0[n];
    bndNodeTag(n)      = bndNodeTag0[n];
    bndNodeNormal(n,0) = bndNodeNormal0[2*n  ];
    bndNodeNormal(n,1) = bndNodeNormal0[2*n+1];
  }
  surfX.allocate(nSurfNode,2);
  for (int n=0; n<nSurfNode; n++){
    surfX(n,0) = surfX0[2*n  ];
    surfX(n,1) = surfX0[2*n+1];
  }
  strandX.allocate(nStrandNode);
  for (int n=0; n<nStrandNode; n++) strandX(n) = strandX0[n];
  pointingVec.allocate(nSurfNode,2);
  for (int n=0; n<nSurfNode; n++){
    pointingVec(n,0) = pointingVec0[2*n  ];
    pointingVec(n,1) = pointingVec0[2*n+1];
  }
  clip.allocate(nSurfNode);
  for (int n=0; n<nSurfNode; n++) clip(n) = clip0[n];


  // if non-corrected scheme is used in the surface direction, convert to
  // a linear mesh
  if (surfOrder <= 2 && meshOrder > 1){
    int meshOrderL=1,nSurfElemL=nSurfElem*meshOrder,nElemEdgeL=meshOrder;
    Array1D<int> surfElemTagL(nSurfElemL);
    Array2D<int>
      surfElemL(nSurfElemL,meshOrderL+1),
      elemEdgeL(nElemEdgeL,2);
    if (meshOrder == 1){
      elemEdgeL(0,0) = 0;
      elemEdgeL(0,1) = 1;
    }
    else if (meshOrder == 2){
      elemEdgeL(0,0) = 0;
      elemEdgeL(0,1) = 2;
      elemEdgeL(1,0) = 2;
      elemEdgeL(1,1) = 1;
    }
    else if (meshOrder == 3){
      elemEdgeL(0,0) = 0;
      elemEdgeL(0,1) = 2;
      elemEdgeL(1,0) = 2;
      elemEdgeL(1,1) = 3;
      elemEdgeL(2,0) = 3;
      elemEdgeL(2,1) = 1;
    }
    else if (meshOrder == 4){
      elemEdgeL(0,0) = 0;
      elemEdgeL(0,1) = 2;
      elemEdgeL(1,0) = 2;
      elemEdgeL(1,1) = 3;
      elemEdgeL(2,0) = 3;
      elemEdgeL(2,1) = 4;
      elemEdgeL(3,0) = 4;
      elemEdgeL(3,1) = 1;
    }

    int k=0;
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<nElemEdgeL; i++){
	surfElemL(k,0)    = surfElem(n,elemEdgeL(i,0));
	surfElemL(k,1)    = surfElem(n,elemEdgeL(i,1));
	surfElemTagL(k++) = surfElemTag(n);
      }
    if (k != nSurfElemL){
      cout << "\n***Problem forming sub-elements in initialize.C***"
	   << endl;
      exit(0);
    }

    meshOrder = meshOrderL;
    nSurfElem = nSurfElemL;
    surfElem.deallocate();
    surfElemTag.deallocate();
    surfElem.allocate(nSurfElem,meshOrder+1);
    surfElemTag.allocate(nSurfElem);
    for (int n=0; n<nSurfElem; n++){
      surfElemTag(n) = surfElemTagL(n);
      for (int i=0; i<meshOrder+1; i++) surfElem(n,i) = surfElemL(n,i);
    }
    surfElemL.deallocate();
    surfElemTagL.deallocate();
    elemEdgeL.deallocate();
  }


  // check to make sure the given mesh order can support the requested
  // solver order of accuracy
  bool problem=false;
  if (meshOrder == 1){
    if (surfOrder > 2) problem = true;
  }
  else if (meshOrder == 2){
    if ((viscous == 0 && surfOrder > 3) ||
	(viscous != 0 && surfOrder > 2)) problem = true;
  }
  else if (meshOrder == 3){
    if (surfOrder > 3) problem = true;
  }
  else if (meshOrder == 4){
    if (surfOrder > 3) problem = true;
  }
  if (problem){
    cout << "\n*** Choose a lower surfOrder for this meshOrder.***"
	 << endl;
    exit(0);
  }

  
  // form all grid connectivity
  connectivity();


  // compute all grid metrics
  metric();


  // allocate and initialize solution variables
  prepare();
}
