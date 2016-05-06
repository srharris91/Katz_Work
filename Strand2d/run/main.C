// File:        main.C
// Release:     $Name:  $
// Revision:    $Revision:  $
// Modified:    $Date:  $
// Description: Driver code to generate and solve on strand meshes

#include "STRANDGEN_defs.h"
#include "STRAND_defs.h"
#include "StrandGlobalMesh.h"
#include "StrandMultiBlockMesh.h"
#include "StrandMultiBlockSolver.h"
#include <stdlib.h>


int main() {
   
   cout << "\n----------------------------------------"
        << "\n         STRANDGEN STRAND DRIVER"
        << "\n            v0.1 July 2010"
        << "\n----------------------------------------"
        << endl;


   // program imputs
   string strand_input_file = "StrandGen2d.input";
   string solver_input_file = "Strand2d.input";


   // create mesh generator instance, read global mesh, and partition
   StrandGlobalMesh::createManager();
   StrandGlobalMesh* sgm = StrandGlobalMesh::getManager();
   sgm->initialize(strand_input_file);


   // invoke the multi-block mesh manager to sprout the mesh
   StrandMultiBlockMesh::createManager();
   StrandMultiBlockMesh* smbm = StrandMultiBlockMesh::getManager();
   smbm->initialize();
   smbm->partition();
   smbm->sprout();
   smbm->print();
   smbm->plot();


   // create multi-block solver, initialize it, and solve for solution
   StrandMultiBlockSolver::createManager();
   StrandMultiBlockSolver* smbs = StrandMultiBlockSolver::getManager();
   smbs->initialize(solver_input_file);
   int restartStep   = smbs->getRestartStep();
   int nSteps        = smbs->getNSteps();
   int nPseudoSteps0 = smbs->getNPseudoSteps0();
   int nPseudoSteps  = smbs->getNPseudoSteps();
   int nPseudoStepsN = nPseudoSteps0;
   int converged     = 0;

   for (int step=restartStep; step<=restartStep+nSteps; step++){
     if (step > 0) nPseudoStepsN = nPseudoSteps;
     for (int pseudoStep=0; pseudoStep<nPseudoStepsN; pseudoStep++){
       smbs->takePseudoStep(step,
			    pseudoStep,
			    converged);
       if (converged == 1 || pseudoStep == nPseudoStepsN-1){
	 smbs->outputStep(step);
	 break;
       }}}
   smbs->finalize();

   
   // free the singleton intances used
   StrandGlobalMesh::freeManager();
   StrandMultiBlockMesh::freeManager();
   StrandMultiBlockSolver::freeManager();
}
