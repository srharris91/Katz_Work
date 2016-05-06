#include "Tri2dFCManager.h"


void Tri2dFCManager::initSolver(string& inputFile)
{
  // for now just have one block solver
  nBlockSolvers = 1;


  // read block solver inputs
  t2dfcbs       = new Tri2dFCBlockSolver[1];
  t2dfcbs[0].input(inputFile);
  nLevels       = t2dfcbs[0].getNLevels();
  iConvFile     = t2dfcbs[0].getIConvFile();
  restartStep   = t2dfcbs[0].getRestartStep();
  nOutput       = t2dfcbs[0].getNOutput();
  nSteps        = t2dfcbs[0].getNSteps();
  nPseudoSteps0 = t2dfcbs[0].getNPseudoSteps0();
  nPseudoSteps  = t2dfcbs[0].getNPseudoSteps();
  nRKStages     = t2dfcbs[0].getNRKStages();
  mgCycle       = t2dfcbs[0].getMgCycle();
  convLimit     = t2dfcbs[0].getConvLimit();
  int* ordersT  = t2dfcbs[0].getOrders();
  orders.allocate(nLevels);
  for (int n=0; n<nLevels; n++) orders(n) = ordersT[n];
  delete [] t2dfcbs;


  // allocate array of block solvers
  t2dfcbs = new Tri2dFCBlockSolver[nBlockSolvers*nLevels];
  cout << "\nSuccessfully allocated " << nBlockSolvers*nLevels
       << " block solvers.\n" << endl;


  // fill in mgLevel and mgMode vectors
  mgMap();


  // for now, no partitioning, just initialize the solver with the global mesh
  for (int level=0; level<nLevels; level++){
    t2dfcbs[level].input(inputFile);
    t2dfcbs[level].initialize(level,
			      orders(level),
			      nTri,
			      nNode,
			      nCompBd,
			      nEdgeBd,
			      tri,
			      &x(0,0),
			      &edgeBd(0,0));
  }

  nq    = t2dfcbs[0].getNq();
  nDofs = t2dfcbs[0].getNDofs();


  // initialize the multigrid transfer operators
  for (int level=1; level<nLevels; level++){
    int     nNodeF   = t2dfcbs[level-1].getNNode();
    int     nNodeBdF = t2dfcbs[level-1].getNNodeBd();
    int*    elemF    = t2dfcbs[level-1].getElem();
    int*    nodeBdF  = t2dfcbs[level-1].getNodeBd();
    double* xF       = t2dfcbs[level-1].getX();
    double* qF       = t2dfcbs[level-1].getQ();
    double* rF       = t2dfcbs[level-1].getR();
    t2dfcbs[level].initRestrict(level-1,
				nNodeF,
				nNodeBdF,
				elemF,
				nodeBdF,
				xF,
				qF,
				rF);
  }

  for (int level=0; level<nLevels-1; level++){
    int*    elemC = t2dfcbs[level+1].getElem();
    double* q0C   = t2dfcbs[level+1].getQ0();
    double* qC    = t2dfcbs[level+1].getQ();
    t2dfcbs[level].initProlong(level+1,
			       elemC,
			       q0C,
			       qC);
  }


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
