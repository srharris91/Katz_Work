/**
 * \mainpage notitle
 * Welcome to the documentation for Tri2dFC2d.
 */

/**
 * \brief
 * Class Tri2dFCBlockSolver holds the data and specifies the operations that 
 * will be carried out on a partitioned block of triangular mesh.  When 
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


#ifndef included_Tri2dFCBlockSolver
#define included_Tri2dFCBlockSolver

#include "TRI2DFC_defs.h"
#include "Tri2dFCSystem.h"


class Tri2dFCBlockSolver
{
 public:
  
  /**
   * \brief
   * Constructor for Tri2dFCBlockSolver.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/Tri2dFCBlockSolver.C Tri2dFCBlockSolver
   */
  Tri2dFCBlockSolver();

  /**
   * \brief
   * Destructor for Tri2dFCBlockSolver
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/Tri2dFCBlockSolver.C ~Tri2dFCBlockSolver
   */
  ~Tri2dFCBlockSolver();

  /**
   * \brief
   * Read and store solver inputs from the input file.
   * \param inputFile name of strand solver input file.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/input.C
   */
  void input(const string& inputFile);

  /**
   * \brief
   * Initializes data structures for the block solver.
   * \param meshFile Meshfile containing unstructured mesh.
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
  void initialize(const int& level,
		  const int& order,
		  const int& nTri,
		  const int& nNode,
		  const int& nCompBd,
		  const int& nEdgeBd,
		  int** tri,
		  const double* x,
		  const int* edgeBd);

  /**
   * \brief
   * Initializes multigrid restriction operator.
   * \param levelF Next fine level number.
   * \param nNodeF Number of nodes on fine level.
   * \param nNodeBdF Number of boundary nodes on fine level.
   * \param elemF Fine level element array.
   * \param xF Fine level coordinates.
   * \param qF Fine level Q array.
   * \param rF Fine level residual array.
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
  void initRestrict(int levelF,
		    int nNodeF,
		    int nNodeBdF,
		    int* elemF,
		    int* nodeBdF,
		    double* xF,
		    double* qF,
		    double* rF);

  /**
   * \brief
   * Initializes multigrid prolongation operator.
   * \param levelC Next coarse level number.
   * \param elemC Coarse level element array.
   * \param q0C Coarse level Q0 array.
   * \param qC Coarse level Q array.
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
  void initProlong(int levelC,
		   int* elemC,
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
   * Deallocates all Tri2dFCBlockSolver data.
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
  const int&    getNLevels();
        int*    getOrders();
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
  const int&    getNNode();
  const int&    getNNodeBd();
  const double& getConvLimit();
        int*    getElem();
        int*    getNodeBd();
        double* getRms();
        double* getX();
        double* getQ();
        double* getQ0();
        double* getR();


 private:

  /**
   * \brief
   * Provides initial values for all Tri2dFCBlockSolver data.
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
   * Form local grid based on global grid.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-05-06
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/localGrid.C
   */
  void localGrid();

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
   * Extract edge data structure for viscous terms.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-05-06
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/edgeExtractVis.C
   */
  void edgeExtractVis();

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
   * Form least squares gradient stencils and coefficients.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2014-03-10
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/leastSquares.C
   */
  void leastSquares();

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
   * Compute nodal volumes.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-05-09
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/volume.C
   */
  void volume();

  /**
   * \brief
   * Compute median dual face areas.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-05-09
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/faceArea.C
   */
  void faceArea();

  /**
   * \brief
   * Compute median dual face areas on each triangular element for
   * viscous discretization.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-05-09
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/faceAreaVis.C
   */
  void faceAreaVis();

  /**
   * \brief
   * Compute FEM gradient coefficients.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-05-09
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/gradSetup.C
   */
  void gradSetup();

  /**
   * \brief
   * Compute quadratic FEM gradient coefficients.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-05-09
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/gradSetupQuadratic.C
   */
  void gradSetupQuadratic();

  /**
   * \brief
   * Compute cubic FEM gradient coefficients.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-05-09
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/gradSetupCubic.C
   */
  void gradSetupCubic();

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
   * Returns solution point locations in the standard equilateral triangle.
   * \param order Order of approximating polynomial.
   * \param spacing Type of point spacing (0=equally spaced, 1=quadrature points).
   * \param rs r and s locations in the standard equilateral triangle.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/solutionPoints.C
   */
  void solutionPoints(const int& order,
		      const int& spacing,
		      double* rs);

  /**
   * \brief
   * Returns solution point locations in the standard equilateral triangle.
   * \param order Order of approximating polynomial.
   * \param triC connectivity in the standard equilateral triangle.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/elementConn.C
   */
  void elementConn(const int& order,
		   int* triC);

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
   * Computes the gradient of a vector.
   * \param npts Number of nodes at which to compute gradients.
   * \param nqp Number of components in the vector.
   * \param p Vector whose gradient is needed.
   * \param px X/Y-components of gradient of p.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-03-26
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/gradient.C
   */
  void gradient(const int& nqp,
		const double* p,
		double* px);

  /**
   * \brief
   * Computes the gradient and Hessian of a vector.
   * \param npts Number of nodes at which to compute Hessian.
   * \param nqp Number of components in the vector.
   * \param p Vector whose Hessian is needed.
   * \param pxx X/Y-components of gradient of p and
   * XX/XY/YY-components of Hessian of p.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-03-26
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/hessian.C
   */
  void hessian(const int& nqp,
	       const double* p,
	       double* pxx);

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


  // Tri2dFC2d solver numerics inputs
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
  static int    nLevels;
  Array1D<int>  orders;
  static int    spacing;
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
  int    order;
  int    nq;
  int    nqa;
  int    ndim;
  int    nTriG;
  int    nNodeG;
  int    nEdgeBdG;
  int    nne;
  int    nee;
  int    nEdgeQ;
  int    nsq;
  int    npsp1;
  int    nElem;
  int    nTri;
  int    nNode;
  int    nNodeBd;
  int    nEdge;
  int    nEdgeBd;
  int    nCompBd;
  int    nqGradQ;
  int    nqaGradQa;
  int    inviscid;
  int    dissipation;
  int    viscous;
  int    source;
  int    sourceMMS;
  int    gradQComputed;
  int    gradQaComputed;
  int    limiterComputed;
  int    nOutputVars;
  int    nOutputScalars;
  int    nOutputVectors;
  double forcex;
  double forcey;
  Tri2dFCSystem* sys;

  int** triG;
  Array2D<double> xG;
  Array2D<int> edgeBdG;

  int     levelC;
  int     levelF;
  int     nNodeF;
  int     nNodeBdF;
  int     nneC;
  int     nneF;
  int*    elemC;
  int*    elemF;
  double* q0C;
  double* qC;
  double* qF;
  double* rF;
  int* nfn;
  int** nfn1;
  double** nfn2;
  Array2D<double> lqCF;
  Array2D<double> lqFC;
  Array2D<int> nfe;
  Array2D<int> nce;

  Array1D<int> iqgrad;
  Array1D<int> iqagrad;
  Array2D<int> elem;
  Array2D<int> tri;
  Array2D<int> edge;
  Array2D<int> edgeE;
  Array2D<int> edgeBd;
  Array1D<int> nodeBd;
  Array2D<int> edgeQ;
  Array1D<int> psp1;
  Array1D<int> psp2;
  Array2D<double> x;
  Array1D<double> v;
  Array2D<double> area;
  Array2D<double> areaBd;
  Array3D<double> areaE;
  Array1D<double> rka;
  Array1D<double> rkb;
  Array1D<double> bdf;
  Array1D<double> dlim;
  Array1D<double> rmsNorm;
  Array1D<double> rms;
  Array2D<double> ln;
  Array4D<double> gxS;
  Array3D<double> nxQ;
  Array4D<double> dxg;
  Array2D<double> gx;
  Array2D<double> gxQ;
  Array2D<double> gxC;
  Array2D<double> gxx;
  Array1D<int> outputVarLength;
  Array1D<string> outputVars;
  Array2D<double> q;
  Array2D<double> qa;
  Array2D<double> q0;
  Array2D<double> fwc;
  Array2D<double> qn;
  Array3D<double> qt;
  Array1D<double> radi;
  Array1D<double> radv;
  Array1D<double> dt;
  Array2D<double> r;
  Array2D<double> d;
  Array2D<double> dn;
  Array2D<double> s;
  Array3D<double> qx;
  Array2D<double> lim;
};
#endif
