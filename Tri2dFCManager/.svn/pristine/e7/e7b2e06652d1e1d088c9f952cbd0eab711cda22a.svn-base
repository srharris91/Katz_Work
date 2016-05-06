// File:        main.C
// Release:     $Name:  $
// Revision:    $Revision:  $
// Modified:    $Date:  $
// Description: Driver code to manage 2d triangle meshes


#include "TRI2DFCMAN_defs.h"
#include "Tri2dFCManager.h"


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
  Tri2dFCManager::createManager();
  Tri2dFCManager* man = Tri2dFCManager::getManager();


  // read manager inputs and mesh file
  man->initialize(inputFile);


  // plot mesh
  man->plot();


  // finalize the manager
  man->finalize();
  Tri2dFCManager::freeManager();
}
