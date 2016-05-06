/**
 * \mainpage notitle
 * Welcome to the documentation for Strand2dFC2d.
 */

/**
 * \brief
 * Class Strand2dFCBlockSolver holds the data and specifies the operations that 
 * will be carried out on a partitioned block of strand mesh.  When 
 * initialized it allocates data and creates communicators between data
 * that needs to be exchanged between blocks.  The other operations
 * define things that need to be done to compute a solution on each block.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2013-01-01
 */


#ifndef included_Strand2dFCBlockSolver
#define included_Strand2dFCBlockSolver

#include "STRAND2DFC_defs.h"
#include "Strand2dFCSystem.h"


class Strand2dFCBlockSolver
{
 public:
  
  /**
   * \brief
   * Constructor for Strand2dFCBlockSolver.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/Strand2dFCBlockSolver.C Strand2dFCBlockSolver
   */
  Strand2dFCBlockSolver();

  /**
   * \brief
   * Destructor for Strand2dFCBlockSolver
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/Strand2dFCBlockSolver.C ~Strand2dFCBlockSolver
   */
  ~Strand2dFCBlockSolver();

  /**
   * \brief
   * Initializes data structures for the block solver.
   * \param inputFile name of strand solver input file.
   * \param level Multigrid level for this block.
   * \param meshOrder Order of the mesh on this block and level.
   * \param nSurfElem Number of surface elements.
   * \param nSurfNode Number of surface nodes.
   * \param nBndNode Number of boundary nodes.
   * \param nStrandNode Number of nodes along a strand.
   * \param nStrandElem Number of elements along a strand.
   * \param nCompBd Number of boundary components.
   * \param surfElem Surface element list.
   * \param bndNode Boundary node list.
   * \param surfX Coordinates of surface nodes.
   * \param strandX 1-d strand coordinates.
   * \param surfElemTag Boundary tag for surface elements.
   * \param bndNodeTag Boundary tag for boundary nodes.
   * \param bndNodeNormal Outward pointing normal at open boundaries
   * \param pointingVec Strand pointing vector.
   * \param clip Clipping index.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-01-01
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/initialize.C
   */
  void initialize(const string& inputFile,
		  const int& level,
		  const int& meshOrder,
		  const int& nSurfElem,
		  const int& nSurfNode,
		  const int& nBndNode,
		  const int& nStrandNode,
		  const int& nCompBd,
		  const int* surfElem,
		  const int* bndNode,
		  const double* surfX,
		  const double* strandX,
		  const int* surfElemTag,
		  const int* bndNodeTag,
		  const double* bndNodeNormal,
		  const double* pointingVec,
		  const int* clip);

  /**
   * \brief
   * Initializes multigrid restriction operator.
   * \param meshOrderF Next fine level mesh order.
   * \param surfElemF Next fine level number of surface elements.
   * \param bndNodeF Next fine level number of boundary nodes.
   * \param qF Next fine level Q array.
   * \param rF Next fine level residual array.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/initRestrict.C
   */
  void initRestrict(const int& meshOrderF,
		    int* surfElemF,
		    int* bndNodeF,
		    double* qF,
		    double* rF);

  /**
   * \brief
   * Initializes multigrid prolongation operator.
   * \param meshOrderC Next coarse level mesh order.
   * \param surfElemC Next coarse level number of surface elements.
   * \param bndNodeC Next coarse level number of boundary nodes.
   * \param q0C Next coarse level Q0 array.
   * \param qC Next coarse level Q array.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/initProlong.C
   */
  void initProlong(const int& meshOrderC,
		   int* surfElemC,
		   int* bndNodeC,
		   double* q0C,
		   double* qC);

  /**
   * \brief
   * Shifts the solution in time one time step for unsteady simulations.
   * \param step The current physical time step number.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-03-26
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/shiftTime.C
   */
  void shiftTime(const int& step);

  /**
   * \brief
   * Computes the full right-hand side discretization.
   * \param step The current physical time step number.
   * \param pseudoStep The current pseudo-time step number.
   * \param stage The current Runge-Kutta stage.
   * \param mode Flag indicating the location withing the multigrid cycle
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-03-26
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/computeRHS.C
   */
  void computeRHS(const int& step,
		  const int& pseudoStep,
		  const int& stage,
		  const int& mode);
  
  /**
   * \brief
   * Solves the linear system and update solution.
   * \param step The current physical time step number.
   * \param pseudoStep The current pseudo-time step number.
   * \param stage The current Runge-Kutta stage.
   * \param mode Flag indicating the location withing the multigrid cycle
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-03-26
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/solve.C
   */
  void solve(const int& step,
	     const int& pseudoStep,
	     const int& stage,
	     const int& mode);

  /**
   * \brief
   * Restricts solution and residual.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-03-26
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/restrict.C
   */
  void restrict();

  /**
   * \brief
   * Solves the linear system and update solution.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-03-26
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/prolong.C
   */
  void prolong();

  /**
   * \brief
   * Output solution, residual, error, and surface data to file for a given
   * time step.
   * \param step Unsteady time step.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-01-01
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/output.C
   */
  void output(const int& nBlocks,
	      const int& step);

  /**
   * \brief
   * Deallocates all Strand2dFCBlockSolver data.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/finalize.C
   */
  void finalize();


  // get methods
  const int&    getIConvFile();
  const int&    getRestartStep();
  const int&    getNOutput();
  const int&    getNSteps();
  const int&    getNPseudoSteps();
  const int&    getNPseudoSteps0();
  const int&    getNRKStages();
  const int&    getMgCycle();
  const int&    getNq();
  const int&    getNDofs();
  const int&    getMeshOrder0();
  const double& getConvLimit();
        int*    getSurfElem0();
        int*    getBndNode();
        int*    getClip();
        double* getRms();
        double* getQ();
        double* getQ0();
        double* getR();


 private:

  /**
   * \brief
   * Provides initial values for all Strand2dFCBlockSolver data.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/dataInit.C
   */
  void dataInit();

  /**
   * \brief
   * Form all grid connectivity.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/connectivity.C
   */
  void connectivity();

  /**
   * \brief
   * Form edge data structure.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-05-06
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/edgeExtract.C
   */
  void edgeExtract();

  /**
   * \brief
   * Reorder nodes with boundary nodes at the end of the list
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-05-06
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/reorderNodes.C
   */
  void reorderNodes();

  /**
   * \brief
   * Form FEM gradient stencils.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-05-06
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/gradStencil.C
   */
  void gradStencil();

  /**
   * \brief
   * Form all grid metrics.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/metric.C
   */
  void metric();

  /**
   * \brief
   * Define 1d finite difference coefficients.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-05-09
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/fdCoeff.C
   */
  void fdCoeff();

  /**
   * \brief
   * Compute mapping related terms.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-05-09
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/mapping.C
   */
  void mapping();

  /**
   * \brief
   * Compute unit normal vectors at each boundary node.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-05-09
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/bNormal.C
   */
  void bNormal();

  /**
   * \brief
   * Setup gradient weight coefficients.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-11-20
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/gradSetup.C
   */
  void gradSetup();

  /**
   * \brief
   * Setup sub-element gradient weight coefficients.
   * \param gx Temp array containing sub-element gradient weights.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-11-20
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/gradSetupSub.C
   */
  void gradSetupSub(Array3D<double>& gx);

  /**
   * \brief
   * Setup full element gradient weight coefficients.
   * \param gx Temp array containing gradient weights.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-11-20
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/gradSetupFull.C
   */
  void gradSetupFull(Array3D<double>& gx);

  /**
   * \brief
   * Allocate and initialize solution variables.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-01-01
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/prepare.C
   */
  void prepare();

  /**
   * \brief
   * Adds unsteady time terms to the RHS residual.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-03-26
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/rhsTime.C
   */
  void rhsTime();

  /**
   * \brief
   * Adds artificial dissipation terms to the RHS residual.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-03-26
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/rhsDissipation.C
   */
  void rhsDissipation();

  /**
   * \brief
   * Computes the gradient of a generic vector, p.
   * \param p Vector whose gradient is to be computed.
   * \param px Gradient of p.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-11-20
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/gradient.C
   */
  void gradient(Array3D<double>& p,
		Array4D<double>& px);

  /**
   * \brief
   * Computes shock capturing limiter.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-03-26
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/limit.C
   */
  void limit();

  /**
   * \brief
   * Adds inviscid terms to the RHS residual.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-03-26
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/inviscid.C
   */
  void rhsInviscid();

  /**
   * \brief
   * Adds viscous terms to the RHS residual.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-03-26
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/rhsViscous.C
   */
  void rhsViscous();

  /**
   * \brief
   * Adds physical source terms to the RHS residual.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-03-26
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/rhsSource.C
   */
  void rhsSource();

  /**
   * \brief
   * Adds MMS source terms to the RHS residual.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-03-26
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/rhsSourceMMS.C
   */
  void rhsSourceMMS();

  /**
   * \brief
   * Adds multigrid source terms to the RHS residual.
   * \param mode Solver mode.
   * \param stage Current Runge-Kutta stage number.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-03-26
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/rhsSourceMG.C
   */
  void rhsSourceMG(const int& mode,
  		   const int& stage);

  /**
   * \brief
   * Modifies residual at boundry nodes to accommodate boundary conditions.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-03-26
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/rhsBoundary.C
   */
  void rhsBoundary();

  /**
   * \brief
   * Compute pseudo-time step.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-03-26
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/pseudoTime.C
   */
  void pseudoTime();

  /**
   * \brief
   * Compute inviscid spectral radius.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-03-29
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/specRadi.C
   */
  void specRadi();

  /**
   * \brief
   * Compute viscous spectral radius.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-03-29
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/specRadv.C
   */
  void specRadv();

  /**
   * \brief
   * Update nodal values.
   * \param step Physical time step.
   * \param stage Runge-Kutta stage for the pseudo-time.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-05-09
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/update.C
   */
  void update(const int& step,
  	      const int& stage);


  // Strand2dFC2d solver numerics inputs
  static int    iPrint;
  static int    iTest;
  static int    iDebug;
  static int    iConvFile;
  static int    iSolnFile;
  static int    iResdFile;
  static int    iErrFile;
  static int    iSurfFile;
  static int    standAlone;
  static int    restartStep;
  static int    nRestart;
  static int    nOutput;
  static int    nSteps;
  static int    nPseudoSteps;
  static int    nPseudoSteps0;
  static int    nLinearSteps;
  static int    nRKStages;
  static int    implicit;
  static int    gradMethod;
  static int    mgCycle;
  static int    limiter;
  static int    timeAcc;
  static double dtUnsteady;
  static double cfl;
  static double vnn;
  static double smooth;
  static double convLimit;
  static double relax;

  // local solver data for this block
  int    ID;
  int    level;
  int    meshOrder;
  int    meshOrder0;
  int    meshOrder0C;
  int    meshOrder0F;
  int    surfOrder;
  int    strandOrder;
  int    nSurfElem;
  int    nSurfElem0;
  int    nSurfNode;
  int    nStrandNode;
  int    nBndNode;
  int    nSurfEdge;
  int    nElemEdge;
  int    nQuadPoint;
  int    npsp1;
  int    npsp1MGS;
  int    npsp1MGR;
  int    npsp1MGC;
  int    nDofs;
  int    nq;
  int    nqa;
  int    ndim;
  int    nCompBd;
  int    nqGradQ;
  int    nqaGradQa;
  int    inviscid;
  int    dissipation;
  int    viscous;
  int    source;
  int    sourceMMS;
  int    nOutputVars;
  int    nOutputScalars;
  int    nOutputVectors;
  int    nDci1;
  int    nDcbEdge1;
  int    nDcb1;
  int    nDci;
  int    nDcbEdge;
  int    nDcb;
  int    nIci;
  int    nIcbNode;
  int    nIcb;
  int    nVci;
  int    nVcbNode;
  int    nVcb;
  int    nVcb2;
  double deltaS;
  double deltaN;
  Strand2dFCSystem* sys;
  int* surfElem0C;
  int* surfElem0F;
  int* bndNodeC;
  int* bndNodeF;
  double* qC;
  double* q0C;
  double* qF;
  double* rF;

  Array2D<int> surfElem;
  Array2D<int> surfElem0;
  Array1D<int> surfElemTag;
  Array2D<int> surfNodeTag;
  Array1D<int> bndNode;
  Array2D<int> bndElem;
  Array1D<int> bndNodeTag;
  Array2D<double> bndNodeNormal;
  Array2D<double> surfX;
  Array1D<double> strandX;
  Array2D<double> pointingVec;
  Array1D<int> clip;
  Array1D<int> iqgrad;
  Array1D<int> iqagrad;
  Array1D<double> rka;
  Array1D<double> rkb;
  Array1D<double> bdf;
  Array1D<double> dlim;
  Array1D<double> rmsNorm;
  Array1D<double> rms;
  Array1D<int> outputVarLength;
  Array1D<string> outputVars;
  Array2D<int> elemEdge;
  Array2D<int> surfEdge;
  Array1D<double> bndSign;
  Array1D<int> psp1;
  Array1D<int> psp2;
  Array1D<int> psp1MGS;
  Array1D<int> psp2MGS;
  Array1D<double> wsp1MGS;
  Array1D<int> psp1MGR;
  Array1D<int> psp2MGR;
  Array1D<double> wsp1MGR;
  Array1D<int> psp1MGC;
  Array1D<int> psp2MGC;
  Array1D<double> wsp1MGC;
  Array2D<double> ls;
  Array2D<double> xn;
  Array2D<double> yn;
  Array3D<double> xs;
  Array3D<double> ys;
  Array2D<double> xsA;
  Array2D<double> ysA;
  Array3D<double> xns;
  Array3D<double> yns;
  Array3D<double> jac;
  Array2D<double> v;
  Array3D<double> sn;
  Array1D<double> wQ;
  Array2D<double> lQ;
  Array2D<double> lsQ;
  Array2D<double> xsQ;
  Array2D<double> ysQ;
  Array2D<double> xnQ;
  Array2D<double> ynQ;
  Array2D<double> jacQ;
  Array3D<double> gxc;
  Array1D<double> dci1;
  Array2D<double> dcb1;
  Array1D<double> dci;
  Array2D<double> dcb;
  Array1D<double> ici;
  Array2D<double> icb;
  Array2D<double> vci;
  Array3D<double> vcb;
  Array2D<int> vciIndex;
  Array3D<int> vcbIndex;
  Array3D<double> q;
  Array3D<double> qa;
  Array3D<double> q0;
  Array3D<double> fwc;
  Array3D<double> qn;
  Array5D<double> qt;
  Array2D<double> radi;
  Array2D<double> radv;
  Array2D<double> dt;
  Array3D<double> r;
  Array3D<double> d;
  Array3D<double> dn;
  Array3D<double> s;
  Array4D<double> qx;
  Array3D<double> lim;
  Array3D<double> strandLim;
};
#endif
