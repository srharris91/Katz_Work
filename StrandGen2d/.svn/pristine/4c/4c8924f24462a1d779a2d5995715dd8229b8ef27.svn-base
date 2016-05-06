// File:        StrandGlobalMesh.h
// Package:     STRANDGEN
//              Parallel Infrastructure for Cartesian & Strand Solvers
// Release:     $Name:  $
// Revision:    $Revision:  $
// Modified:    $Date:  $
// Description: Data structure holding global mesh; partitions mesh and
//              creates StrandBlocks


#ifndef included_StrandGlobalMesh
#define included_StrandGlobalMesh

#include "STRANDGEN_defs.h"
#include "StrandBlock.h"


class StrandGlobalMesh
{
public:

  // @brief Create singleton instance.
  static void createManager();

  // @brief Return a pointer to the singleton instance.
  static StrandGlobalMesh* getManager();

  // @brief Deallocate the manager instance.
  static void freeManager();

  // initialize the global mesh (read input and mesh files, allocate space for
  // the global mesh and store it)
  void initialize(string& input_file);
 
  // get methods:
  const bool& getInitialized();
  const int& getNPartitions();
  const int& getIwrite();
  const int& getStrandDist();
  const int& getNPtsPerStrand();
  const int& getNFringe();
  const int& getTrimType();
  const int& getSmax();
  const double& getAngle();
  const double& getStretchRatio();
  const double& getWallSpacing();
  const double& getStrandLength();
  const double& getDeltaSmooth();

  const int& getSurfaceOnly();
  const int& getNSurfNodes();
  const int& getNSurfFaces();
  const int& getNBndEdges();
  const int& getNSurfPatches();
  const int& getNEdgePatches();
  const int& getNSharp();
  const int* getSurfFaces();
  const int* getBndEdges();
  const int* getNodeClip();
  const int* getFaceClip();
  const int* getFaceTag();
  const int* getBndTag();
  const int* getSFlag();
  const double* getXSurf();
  const double* getPointingVec();
  const double* getXStrand();
  const double* getBndNormal();

  // set methods:
  void setNPartitions(int);


protected:   

  // @brief protected constructor,
  // protected to ensure only one instance may be created.
  StrandGlobalMesh();

  // @brief protected destructor.
  //De-allocates any objects resident to this one before it is destroyed.
  virtual ~StrandGlobalMesh();
   
  
private:

  static StrandGlobalMesh* sgm_manager_instance;

  static bool    d_initialized;
  static int     d_nPartitions;
  static int     d_iwrite;
  static int     d_strandDist;
  static int     d_nPtsPerStrand;
  static int     d_nFringe;
  static int     d_trimType;
  static int	 d_smax;
  static double  d_angle;
  static double  d_stretchRatio;
  static double  d_wallSpacing;
  static double  d_strandLength;
  static double  d_deltaSmooth;

  static int     d_surfaceOnly;
  static int     d_nMesh;
  static int     d_nSurfNodes;
  static int     d_nSurfFaces;
  static int     d_nBndEdges;
  static int     d_nSurfPatches;
  static int     d_nEdgePatches;
  static int	 d_nSharp;
  static int*    d_surfFaces;
  static int*    d_bndEdges;
  static int*    d_nodeClip;
  static int*    d_faceClip;
  static int*    d_faceTag;
  static int*    d_bndTag;
  static int*	 d_sFlag;
  static double* d_xSurf;
  static double* d_pointingVec;
  static double* d_xStrand;
  static double* d_bndNormal;
};
#endif  // included_StrandGlobalMesh
