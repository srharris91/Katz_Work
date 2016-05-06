/**
 * \mainpage notitle
 * Welcome to the documentation for the Strand2dFC Manager.
 */

/**
 * \brief
 * Class Strand2dFCManager is the global mesh and solver manager for the
 * Strand2dFC flux correction code on strands.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 03-15-2013
 */


#ifndef included_Strand2dFCManager
#define included_Strand2dFCManager

#include "STRAND2DFCMAN_defs.h"
#include "Strand2dFCBlockMesh.h"
#include "Strand2dFCBlockSolver.h"


class Strand2dFCManager
{
 public:

  /**
   * \brief
   * Creates a singleton instance object for class Strand2dFCManager.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-15-2013
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Strand2dFCManager.C createManager
   */
  static void createManager();


  /**
   * \brief
   * Returns a pointer to the singleton instance object for
   * class Strand2dFCManager.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-15-2013
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Strand2dFCManager.C getManager
   */
  static Strand2dFCManager* getManager();


  /**
   * \brief
   * Frees the pointer to the singleton instance object for
   * class Strand2dFCManager.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-15-2013
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Strand2dFCManager.C freeManager
   */
  static void freeManager();


  /**
   * \brief
   * Initializes the mesh manager.
   * \param inputFile Name of input file.
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
   * Creates an output file (Paraview format) of the mesh.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-15-2013
   * \par Further Documentation:
   * \par Source Code:
   */
  void plot();

  /**
   * \brief
   * Initializes the block solvers.
   * \param inputFile Name of input file.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-18-2013
   * \par Further Documentation:
   * \par Source Code:
   * \include src/initializeSolver.C
   */
  void initializeSolver(string& inputFile);


  /**
   * \brief
   * Take pseudo-time step.
   * \param step Current physical time step.
   * \param pseudoStep Current pseudo time step.
   * \param converged Flag to determine if converged.
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
   * \param step Current physical time step.
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
   * Deallocates memory and finalizes the block solvers.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-15-2013
   * \par Further Documentation:
   * \par Source Code:
   * \include src/finalizeSolver.C
   */
  void finalizeSolver();


  // get methods
  const int&    getRestartStep();
  const int&    getNOutput();
  const int&    getNSteps();
  const int&    getNPseudoSteps();
  const int&    getNPseudoSteps0();


 protected:   

  /**
   * \brief
   * Protected constructor for class Strand2dFCManager.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-15-2013
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Strand2dFCManager.C Strand2dFCManager
   */
  Strand2dFCManager();

  /**
   * \brief
   * Protected destructor for class Strand2dFCManager.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 03-15-2013
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Strand2dFCManager.C ~Strand2dFCManager
   */
  virtual ~Strand2dFCManager();


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
   * Report convergence history.
   * \param step Current physical time step.
   * \param pseudoStep Current pseudo time step.
   * \param converged Flag to determine if converged.
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


  static Strand2dFCManager* Strand2dfc_manager_instance;
  static int nLevels;
  static int iplotmesh;
  static int nBlockSolvers;
  static int nBlockMeshes;
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
  ofstream cfile;
  Array2D<Strand2dFCBlockSolver> blockSolver;
  Array2D<Strand2dFCBlockMesh> blockMesh;
  Array1D<int> mgLevel;
  Array1D<int> mgMode;
};
#endif
