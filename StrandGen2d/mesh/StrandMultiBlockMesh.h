// File:        StrandMultiBlockMesh.h
// Package:     STRANDGEN
//              Parallel Infrastructure for Cartesian & Strand Solvers
// Release:     $Name:  $
// Revision:    $Revision:  $
// Modified:    $Date:  $
// Description: Manages strand mesh generation across all blocks

/*!
 * StrandMultiBlockMesh manages the global strand mesh generation by
 * invoking StrandBlock member functions in loops over blocks.
 *
 * The main functions in this class are:
 *    - @b    extract()
 *      Extracts the surface mesh edges.
 *
 *
 * The class is a SINGLETON object, meaning there can be only one instance
 * in the application. This instance can be accessed from any other 
 * object, at any time, using the following convention:
 *
 *    \verbatim
 *    StrandMultiBlockMesh::createManager();    // must be called once
 *    StrandMultiBlockMesh* smbm = StrandMultiBlockMesh::getManager();
 *    smbm->extract(); 
 *    \endverbatim
 */

#ifndef included_StrandMultiBlockMesh
#define included_StrandMultiBlockMesh

#include "STRANDGEN_defs.h"
#include "StrandGlobalMesh.h"
#include "StrandBlock.h"
#include <stdlib.h>

class StrandMultiBlockMesh
{
public:

  // @brief Create singleton instance.
  static void createManager();

  // @brief Return a pointer to the singleton instance.
  static StrandMultiBlockMesh* getManager();

  // @brief Deallocate the manager instance.
  static void freeManager();

  // @brief Initialize the mesh manager
  void initialize();

  // @brief Partition global mesh in to tri blocks
  void partition();

  // @brief Sprout the volume mesh
  void sprout();

  // @brief Print information on blocks
  void print();

  // @brief Create a plot of the mesh
  void plot();

  // get methods
  const int&   getNStrandBlocks();
  const int&   getNSurfFaces();
  const int&   getNSurfNodes();
  const int&   getNBndEdges();
  const int&   getNPrtEdges();
  const int*   getGlobalFaceMap();
  const int*   getGlobalNodeMap();
  const int*   getGlobalEdgeMap();
  const int*   getPrtEdges();
  StrandBlock* getStrandBlock(int&);
  int* getNFaces();
  int* getNNodes();
  int* getNBedges();
  int* getNPstr();
  int* getNFringe();
  int** getFace();
  int** getBEdge();
  double** getXSrf();
  double** getXStr();
  double** getPV();


 protected:

  // @brief protected  constructor.
  StrandMultiBlockMesh();

  // @brief protected destructor.
  virtual ~StrandMultiBlockMesh();
   

private:

  // singleton instance 
  static StrandMultiBlockMesh* smbm_manager_instance;

  static int  d_nStrandBlocks;
  static int  d_nSurfFaces;
  static int  d_nSurfNodes;
  static int  d_nBndEdges;
  static int  d_nPrtEdges;
  static int* d_globalFaceMap;
  static int* d_globalNodeMap;
  static int* d_globalEdgeMap;
  static int* d_prtEdges;
  static StrandBlock* d_strandBlocks;
  static int* d_nFaces;
  static int* d_nNodes;
  static int* d_nBedges;
  static int* d_nPstr;
  static int* d_nFringe;
  static int** d_face;
  static int** d_bEdge;
  static double** d_xSrf;
  static double** d_xStr;
  static double** d_pV;

  // status logicals
  static bool d_initialized;
  static bool d_partitioned;
  static bool d_sprouted;
};
#endif
