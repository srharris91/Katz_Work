// File:        mainSolver.C
// Release:     $Name:  $
// Revision:    $Revision:  $
// Modified:    $Date:  $
// Description: Driver code to solve on 2d strand meshes


#include "STRAND2DFCMAN_defs.h"
#include "Strand2dFCManager.h"


int main(int argc,
	 char *argv[]){

  cout << "\n----------------------------------------"
       << "\n          SOLVER TEST DRIVER"
       << "\n            v0.1 March 2013"
       << "\n----------------------------------------"
       << endl;
  
  if (argc != 2){
    cout << "\nSpecify name of input file in double quotes." << endl;
    exit(0);
  }

  string inputFile = argv[1];

  cout << "\n"
       << "\nRunning executable " << argv[0]
       << "\nUsing input file " << inputFile
       << "\n"
       << endl;


  // create and initialize global manager
  Strand2dFCManager::createManager();
  Strand2dFCManager* man = Strand2dFCManager::getManager();


  // initialize manager and solver
  man->initialize(inputFile);
  man->initializeSolver(inputFile);


  // obtain solver parameters
  int  restartStep   = man->getRestartStep();
  int  nOutput       = man->getNOutput();
  int  nSteps        = man->getNSteps();
  int  nPseudoSteps0 = man->getNPseudoSteps0();
  int  nPseudoSteps  = man->getNPseudoSteps();
  int  nPseudoStepsN = nPseudoSteps0;
  bool converged     = false;


  // solve
  for (int step=restartStep; step<=restartStep+nSteps; step++){
    if (step > 0) nPseudoStepsN = nPseudoSteps;
    for (int pseudoStep=0; pseudoStep<nPseudoStepsN; pseudoStep++){
      man->takePseudoStep(step,
			  pseudoStep,
			  converged);
      if (converged) break;
    }
    if (step % nOutput == 0) man->output(step);
  }


  // finalize the manager and solver
  man->finalizeSolver();
  man->finalize();
  Strand2dFCManager::freeManager();
}
