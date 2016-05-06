/**
 * \brief
 * Class Strand2dFCBlockMesh holds the data and specifies the operations that 
 * will be carried out on a partitioned block of strand mesh.  When 
 * initialized it allocates all mesh data and generates pointing vectors
 * and strand node distributions.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2013-01-01
 */


#ifndef included_Strand2dFCBlockManager
#define included_Strand2dFCBlockManager

#include "STRAND2DFCMAN_defs.h"
#include "solutionPoints1D.h"
#include "lagrangePoly1D.h"


class Strand2dFCBlockMesh
{
 public:
  
  /**
   * \brief
   * Constructor for Strand2dFCBlockMesh.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Strand2dFCBlockMesh.C Strand2dFCBlockMesh
   */
  Strand2dFCBlockMesh();
  
  /**
   * \brief
   * Destructor for Strand2dFCBlockMesh
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Strand2dFCBlockMesh.C ~Strand2dFCBlockMesh
   */
  ~Strand2dFCBlockMesh();
  
  /**
   * \brief
   * Initializes data structures for the block mesh.
   * \param level Multigrid level for this block.
   * \param meshOrder Order of the mesh on this block and level.
   * \param nSurfElemG Number of surface elements in global mesh.
   * \param nSurfNodeG Number of surface nodes in global mesh.
   * \param nBndNodeG Number of boundary nodes in global mesh.
   * \param nStrandNodeG Number of nodes along a strand in global mesh.
   * \param nCompBdG Number of boundary components in global mesh.
   * \param surfElemG Surface element list in global mesh.
   * \param bndNodeG Boundary node list in global mesh.
   * \param surfXG Coordinates of surface nodes in global mesh.
   * \param strandXG Coordinates of strand nodes in global mesh.
   * \param surfElemTagG Boundary tag for surface elements in global mesh.
   * \param bndNodeTagG Boundary tag for boundary nodes in global mesh.
   * \param bndNodeNormalG Outward pointing normal at open boundaries in
   *                       global mesh.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-11-22
   * \par Further Documentation:
   * \par Source Code:
   * \include src/initializeMesh.C
   */
  void initialize(const int& level,
		  const int& meshOrder,
		  const int& nSurfElemG,
		  const int& nSurfNodeG,
		  const int& nBndNodeG,
		  const int& nStrandNodeG,
		  const int& nCompBdG,
		  int** surfElemG,
		  const Array1D<int>& bndNodeG,
		  const Array2D<double>& surfXG,
		  const Array1D<double>& strandXG,
		  const Array1D<int>& surfElemTagG,
		  const Array1D<int>& bndNodeTagG,
		  const Array2D<double>& bndNodeNormalG);
  
  /**
   * \brief
   * Creates an output file (Paraview format) of the mesh.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-11-22
   * \par Further Documentation:
   * \par Source Code:
   * \include src/plotMesh.C
   */
  void plot();

  /**
   * \brief
   * Finalizes the block mesh.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-11-22
   * \par Further Documentation:
   * \par Source Code:
   * \include src/finalize.C
   */
  void finalize();
  
  
  // get methods
  const int&    getMeshOrder();
  const int&    getNSurfElem();
  const int&    getNSurfNode();
  const int&    getNBndNode();
  const int&    getNStrandNode();
  const int&    getNCompBd();
        int*    getSurfElem();
        int*    getBndNode();
        double* getSurfX();
        double* getStrandX();
        int*    getSurfElemTag();
        int*    getBndNodeTag();
        double* getBndNodeNormal();
        double* getPointingVec();
        int*    getClip();

  
 private:
  
  // local mesh data for this block
  int level;
  int meshOrder;
  int nSurfElem;
  int nSurfNode;
  int nBndNode;
  int nStrandNode;
  int nCompBd;
  Array2D<int> surfElem;
  Array2D<double> surfX;
  Array1D<int> bndNode;
  Array1D<int> surfElemTag;
  Array1D<int> bndNodeTag;
  Array2D<double> bndNodeNormal;
  Array2D<double> pointingVec;
  Array1D<int> clip;
  Array1D<double> strandX;
};
#endif
