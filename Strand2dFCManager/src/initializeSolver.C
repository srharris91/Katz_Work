#include "Strand2dFCManager.h"


void Strand2dFCManager::initializeSolver(string& inputFile)
{
  // allocate array of block solvers
  blockSolver.allocate(nBlockSolvers,nLevels);


  // for now, no partitioning, just initialize the solvers with the
  // global mesh data
  for (int level=0; level<nLevels; level++){
    int     meshOrder=blockMesh(0,level).getMeshOrder();
    int     nSurfElem=blockMesh(0,level).getNSurfElem();
    int     nSurfNode=blockMesh(0,level).getNSurfNode();
    int     nBndNode=blockMesh(0,level).getNBndNode();
    int     nStrandNode=blockMesh(0,level).getNStrandNode();
    int     nCompBd=blockMesh(0,level).getNCompBd();
    int*    surfElem=blockMesh(0,level).getSurfElem();
    int*    bndNode=blockMesh(0,level).getBndNode();
    double* surfX=blockMesh(0,level).getSurfX();
    double* strandX=blockMesh(0,level).getStrandX();
    int*    surfElemTag=blockMesh(0,level).getSurfElemTag();
    int*    bndNodeTag=blockMesh(0,level).getBndNodeTag();
    double* bndNodeNormal=blockMesh(0,level).getBndNodeNormal();
    double* pointingVec=blockMesh(0,level).getPointingVec();
    int*    clip=blockMesh(0,level).getClip();
    blockSolver(0,level).initialize(inputFile,
				    level,
				    meshOrder,
				    nSurfElem,
				    nSurfNode,
				    nBndNode,
				    nStrandNode,
				    nCompBd,
				    surfElem,
				    bndNode,
				    surfX,
				    strandX,
				    surfElemTag,
				    bndNodeTag,
				    bndNodeNormal,
				    pointingVec,
				    clip);
  }


  // initialize the multigrid transfer operators
  for (int level=1; level<nLevels; level++){
    int     meshOrder0F=blockSolver(0,level-1).getMeshOrder0();
    int*    surfElem0F=blockSolver(0,level-1).getSurfElem0();
    int*    bndNodeF=blockSolver(0,level-1).getBndNode();
    double* qF=blockSolver(0,level-1).getQ();
    double* rF=blockSolver(0,level-1).getR();
    blockSolver(0,level).initRestrict(meshOrder0F,
				      surfElem0F,
				      bndNodeF,
				      qF,
				      rF);
  }

  for (int level=0; level<nLevels-1; level++){
    int     meshOrder0C=blockSolver(0,level+1).getMeshOrder0();
    int*    surfElem0C=blockSolver(0,level+1).getSurfElem0();
    int*    bndNodeC=blockSolver(0,level+1).getBndNode();
    double* q0C=blockSolver(0,level+1).getQ0();
    double* qC=blockSolver(0,level+1).getQ();
    blockSolver(0,level).initProlong(meshOrder0C,
				     surfElem0C,
				     bndNodeC,
				     q0C,
				     qC);
  }


  // get solver parameters from the block solver
  nq            = blockSolver(0,0).getNq();
  nDofs         = blockSolver(0,0).getNDofs();
  iConvFile     = blockSolver(0,0).getIConvFile();
  restartStep   = blockSolver(0,0).getRestartStep();
  nOutput       = blockSolver(0,0).getNOutput();
  nSteps        = blockSolver(0,0).getNSteps();
  nPseudoSteps0 = blockSolver(0,0).getNPseudoSteps0();
  nPseudoSteps  = blockSolver(0,0).getNPseudoSteps();
  nRKStages     = blockSolver(0,0).getNRKStages();
  mgCycle       = blockSolver(0,0).getMgCycle();
  convLimit     = blockSolver(0,0).getConvLimit();


  // fill in mgLevel and mgMode vectors
  mgMap();


  // initialize on-screen and convergence file output
  if (iConvFile != 0){
    cfile.open("convergence.dat");
    cfile.setf(ios::scientific);
  }
  time0 = clock();
  cout << "\n";
  cout << "step        pseudoStep  time(s)     RMS residuals" << endl;
  cout << "----------  ----------  ----------  -------------" << endl;
}
