// File:        main.C
// Release:     $Name:  $
// Revision:    $Revision:  $
// Modified:    $Date:  $
// Description: Driver code to manage 2d strand meshes


#include "STRAND2DFCMAN_defs.h"
#include "Strand2dFCManager.h"


int main(int argc,
	 char *argv[]){

  cout << "\n----------------------------------------"
       << "\n          MANAGER TEST DRIVER"
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


  // read manager inputs and mesh file
  man->initialize(inputFile);


  // plot mesh
  man->plot();


  // finalize the manager
  man->finalize();
  Strand2dFCManager::freeManager();
}
