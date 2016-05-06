// File:        main.C
// Release:     $Name:  $
// Revision:    $Revision:  $
// Modified:    $Date:  $
// Description: Driver code to generate strand meshes

#include "STRANDGEN_defs.h"
#include "StrandGlobalMesh.h"
#include "StrandMultiBlockMesh.h"
#include <stdlib.h>

int main( int argc, char *argv[] ) {
   
   cout << "\n----------------------------------------"
        << "\n         STRANDGEN STRAND DRIVER"
        << "\n            v0.1 July 2010"
        << "\n----------------------------------------"
        << endl;

   // program imputs
   string mesh_input_file = "StrandGen2d.input";


   // create mesh generator instance, read global mesh, and partition
   StrandGlobalMesh::createManager();
   StrandGlobalMesh* sgm = StrandGlobalMesh::getManager();
   sgm->initialize(mesh_input_file);


   // invoke the multi-block mesh manager to sprout the mesh
   StrandMultiBlockMesh::createManager();
   StrandMultiBlockMesh* smbm = StrandMultiBlockMesh::getManager();
   smbm->initialize();
   smbm->partition();
   smbm->sprout();
   smbm->print();
   smbm->plot();


   // free the singleton intances used
   StrandGlobalMesh::freeManager();
   StrandMultiBlockMesh::freeManager();
}
