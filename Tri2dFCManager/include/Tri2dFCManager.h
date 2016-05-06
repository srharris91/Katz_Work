/**
 * \mainpage notitle
 * Welcome to the documentation for the Tri2dFC Manager.
 */

/**
 * \brief
 * Class Tri2dFCManager is the global mesh and solver manager for the
 * Tri2dFC flux correction code on triangles.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 03-15-2013
 */


#ifndef included_Tri2dFCManager
#define included_Tri2dFCManager

#include "TRI2DFCMAN_defs.h"
#include "Tri2dFCBlockSolver.h"


class Tri2dFCManager
{
 public:

  /**
   * \brief
   * Creates a singleton instance object for class Tri2dFCManager.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-15-2013
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Tri2dFCManager.C createManager
   */
  static void createManager();


  /**
   * \brief
   * Returns a pointer to the singleton instance object for
   * class Tri2dFCManager.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-15-2013
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Tri2dFCManager.C getManager
   */
  static Tri2dFCManager* getManager();


  /**
   * \brief
   * Frees the pointer to the singleton instance object for
   * class Tri2dFCManager.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-15-2013
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Tri2dFCManager.C freeManager
   */
  static void freeManager();


  /**
   * \brief
   * Initializes the manager.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-15-2013
   * \par Further Documentation:
   * \par Source Code:
   * \include src/initialize.C
   */
  void initialize(string& inputFile);


  /**
   * \brief
   * Initializes the block solvers.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-18-2013
   * \par Further Documentation:
   * \par Source Code:
   * \include src/initSolver.C
   */
  void initSolver(string& inputFile);


  /**
   * \brief
   * Take pseudo-time step.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-26-2013
   * \par Further Documentation:
   * \par Source Code:
   * \include src/takePseudoStep.C
   */
  void takePseudoStep(const int& step,
		      const int& pseudoStep,
		      bool& converged);


  /**
   * \brief
   * Creates output files (Paraview format) of the solution.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-24-2013
   * \par Further Documentation:
   * \par Source Code:
   * \include src/output.C
   */
  void output(const int& step);


  /**
   * \brief
   * Deallocates memory and finalizes the manager.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-15-2013
   * \par Further Documentation:
   * \par Source Code:
   * \include src/finalize.C
   */
  void finalize();


  /**
   * \brief
   * Creates an output file (Paraview format) of the mesh.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-15-2013
   * \par Further Documentation:
   * \par Source Code:
   * \include src/plot.C
   */
  void plot();

  // get methods
  const int&    getRestartStep();
  const int&    getNOutput();
  const int&    getNSteps();
  const int&    getNPseudoSteps();
  const int&    getNPseudoSteps0();


 protected:   

  /**
   * \brief
   * Protected constructor for class Tri2dFCManager.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-15-2013
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Tri2dFCManager.C Tri2dFCManager
   */
  Tri2dFCManager();

  /**
   * \brief
   * Protected destructor for class Tri2dFCManager.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-15-2013
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Tri2dFCManager.C ~Tri2dFCManager
   */
  virtual ~Tri2dFCManager();


 private:

  /**
   * \brief
   * Determine multigrid topology.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-28-2013
   * \par Further Documentation:
   * \par Source Code:
   * \include src/mgMap.C
   */
  void mgMap();

  /**
   * \brief
   * Performs multigrid restriction.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-28-2013
   * \par Further Documentation:
   * \par Source Code:
   * \include src/restrict.C
   */
  void restrict(const int& level);

  /**
   * \brief
   * Performs multigrid prolongation.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-28-2013
   * \par Further Documentation:
   * \par Source Code:
   * \include src/prolong.C
   */
  void prolong(const int& level);

  /**
   * \brief
   * Report convergence history.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-28-2013
   * \par Further Documentation:
   * \par Source Code:
   * \include src/convHistory.C
   */
  void convHistory(const int& step,
		   const int& pseudoStep,
		   bool& converged);


  static Tri2dFCManager* tri2dfc_manager_instance;
  static int iplotmesh;
  static int nTri;
  static int nNode;
  static int nCompBd;
  static int nEdgeBd;
  static int nBlockSolvers;
  static int nLevels;
  static int iConvFile;
  static int restartStep;
  static int nOutput;
  static int nSteps;
  static int nPseudoSteps0;
  static int nPseudoSteps;
  static int nRKStages;
  static int mgCycle;
  static int nq;
  static int nDofs;
  static double convLimit;
  static clock_t time0;
  static clock_t timeS0;
  static clock_t time;
  static Tri2dFCBlockSolver* t2dfcbs;
  ofstream cfile;
  static int** tri;
  Array2D<double> x;
  Array2D<int> edgeBd;
  Array1D<int> orders;
  Array1D<int> mgLevel;
  Array1D<int> mgMode;
};
#endif
