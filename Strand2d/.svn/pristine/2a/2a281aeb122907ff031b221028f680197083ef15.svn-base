/**
 * \mainpage notitle
 * Welcome to the documentation for Strand2d.
 */

/**
 * \brief
 * Class StrandBlockSolver holds the data and specifies the operations that 
 * will be carried out on a partitioned block of the strand mesh.  When 
 * initialized it allocates data and creates communicators between data
 * that needs to be exchanged between blocks.  The other operations
 * define things that need to be done to compute a solution on each block.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-08-30
 */


#ifndef included_StrandBlockSolver
#define included_StrandBlockSolver

#include "STRAND_defs.h"
#include "StrandBlock.h"
#include "StrandSystem.h"


class StrandBlockSolver
{
 public:
  
  /**
   * \brief
   * Constructor for StrandBlockSolver.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/StrandBlockSolver.C StrandBlockSolver
   */
  StrandBlockSolver();

  /**
   * \brief
   * Copy constructor for StrandBlockSolver.
   * \param sbs Reference to the StrandBlockSolver for the copy constructor.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/StrandBlockSolver.C StrandBlockSolver
   */
  StrandBlockSolver(const StrandBlockSolver& sbs);

  /**
   * \brief
   * Destructor for StrandBlockSolver
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/StrandBlockSolver.C ~StrandBlockSolver
   */
  ~StrandBlockSolver();

  /**
   * \brief
   * Read and store global inputs from the Strand3d input deck.
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
   * Retrieves data from the StrandBlock and copies it into StrandBlockSolver
   * data.
   * \param block The StrandBlock corresponding to this StrandBlockSolver.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/setStrandBlockData.C
   */
  void setStrandBlockData(StrandBlock& strandBlock);

  /**
   * \brief
   * Allocates all solver data and performs grid and flow initialization.
   * \param mglevel The current multigrid level.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/initialize.C
   */
  void initialize(const int& mglevel);
  
  /**
   * \brief
   * Creates coarse multigrid level data given a parent StrandBlockSolver.
   * \param mglevel The current multigrid level.
   * \param parentSbs The parent StrandBlockSolver from which the current
   * block is coarsened for multigrid.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/coarsen.C
   */
  void coarsen(const int& mglevel,
	       StrandBlockSolver& parentSbs);
  
  /**
   * \brief
   * Shifts the solution in time one time step for unsteady simulations.
   * \param step The current physical time step number.
   * \param mglevel The current multigrid level.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/shiftTime.C
   */
  void shiftTime(const int& step,
		 const int& mglevel);
  
  /**
   * \brief
   * Outputs the solution on the given StrandBlockSolver.
   * \param step The current physical time step number.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/output.C
   */
  void output(const int& step);
  
  /**
   * \brief
   * Compute gradients of q at interior and boundary cells (not fringe).
   * \param mglevel The current multigrid level.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/gradQ.C gradQ
   */
  void gradQ(const int& mglevel);

  /**
   * \brief
   * Compute gradients of qa at interior and boundary cells (not fringe).
   * \param mglevel The current multigrid level.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/gradQa.C gradQa
   */
  void gradQa(const int& mglevel);

  /**
   * \brief
   * Computes the full right-hand side discretization.
   * \param step The current physical time step number.
   * \param pseudoStep The current pseudo-time step number.
   * \param pseudoStep The current linear-time step number.
   * \param sweep The current GS sweep.
   * \param mglevel The current multigrid level.
   * \param mode Flag indicating the location withing the multigrid cycle
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/computeRHS.C
   */
  void computeRHS(const int& step,
		  const int& pseudoStep,
		  const int& linearStep,
		  const int& sweep,
		  const int& mglevel,
		  const int& mode);
  
  /**
   * \brief
   * Computes the left-hand side linearization terms.
   * \param step The current physical time step number.
   * \param pseudoStep The current pseudo-time step number.
   * \param mglevel The current multigrid level.
   * \param mode Flag indicating the location withing the multigrid cycle
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/computeLHS.C
   */
  void computeLHS(const int& step,
		  const int& pseudoStep,
		  const int& mglevel,
		  const int& mode);
  
  /**
   * \brief
   * Solves the linear system and update solution.
   * \param step The current physical time step number.
   * \param pseudoStep The current pseudo-time step number.
   * \param linearStep The current linear-time step number.
   * \param mglevel The current multigrid level.
   * \param mode Flag indicating the location withing the multigrid cycle
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/solve.C
   */
  void solve(const int& step,
	     const int& pseudoStep,
	     const int& linearStep,
	     const int& mglevel,
	     const int& mode);
  
  /**
   * \brief
   * Restricts q and r from a fine MG level to a coarse MG level.
   * \param mglevel The current multigrid level.
   * \param mode Flag indicating the location withing the multigrid cycle
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/restrict.C
   */
  void restrict(const int& mglevel,
		const int& mode);
  
  /**
   * \brief
   * Prolongs q corrections from a coarse MG level to a fine MG level.
   * \param mglevel The current multigrid level.
   * \param mode Flag indicating the location withing the multigrid cycle.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * <a href="../../notes/prolong.pdf" target="_self"><b>prolong.pdf</b></a>
   * \par Source Code:
   * \include src/Numerics/prolong.C
   */
  void prolong(const int& mglevel,
	       const int& mode);
  
  /**
   * \brief
   * Deallocates all StrandBlockSolver data.
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
  void finalize(const int& mglevel);
  

  // get methods
  const int&    getIConvFile();
  const int&    getRestartStep();
  const int&    getNSteps();
  const int&    getNPseudoSteps();
  const int&    getNPseudoSteps0();
  const int&    getNLinearSteps();
  const int&    getNLevels();
  const int&    getMgCycle();
  const int&    getNq();
  const int&    getNFaces();
  const int&    getNPstr();
  const int&    getNBedges();
  const int&    getNDofs();
  const double& getConvLimit();
        int*    getFClip();
  Array1D<int>* getFClipArray();
        double* getQ();
  Array3D<double>* getQArray();
        double* getXc();
        double* getQx();
	double* getV();
  Array2D<double>* getVArray();
        double* getRms();


 private:

  /**
   * \brief
   * Provides initial values for all StrandBlockSolver data.
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
   * Orders cells for GS procedure.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/order.C
   */
  void order();

  /**
   * \brief
   * Form implicit cells surrounding cells connectivity.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-08-30
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/formCsc.C
   */
  void formCsc();

  /**
   * \brief
   * perturb nodal locations for verification testing.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-11-12
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/PerturbNodes.C
   */
  void perturbNodes();

  /**
   * \brief
   * Computes least squares nodal interpolation coefficient via volume
   * weighting.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-20
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/lspVol.C
   */
  void lspVol();

  /**
   * \brief
   * Computes least squares nodal interpolation coefficient via least squares.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-21
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/lspLS.C
   */
  void lspLS();

  /**
   * \brief
   * Computes least squares nodal interpolation coefficient via a mapped
   * least squares procedure.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-21
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/lspMap.C
   */
  void lspMap();

  /**
   * \brief
   * Computes nodal values of q from cell-center values using lsp arrays.
   * \param mglevel Multigrid level.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-08
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/nodalQ.C nodalQ
   */
  void nodalQ(const int& mglevel);

  /**
   * \brief
   * Computes nodal values of qa from cell-center values using lsp arrays.
   * \param mglevel Multigrid level.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-08
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/nodalQa.C nodalQa
   */
  void nodalQa(const int& mglevel);

  /**
   * \brief
   * Computes coarse metrics (volumes and face areas) via agglomeration.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-17
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/coarseMetrics.C
   */
  void coarseMetrics();

  /**
   * \brief
   * Adds physical time derivative terms to the RHS on fine MG levels.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-18
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/rhsTimeFine.C rhsTimeFine
   */
  void rhsTimeFine();

  /**
   * \brief
   * Adds physical time derivative terms to the RHS on coarse MG levels.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-18
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/rhsTimeCoarse.C rhsTimeCoarse
   */
  void rhsTimeCoarse();

  /**
   * \brief
   * Adds dissipation and inviscid terms to the RHS on fine MG levels.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-18
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/rhsDissipationInviscidFine.C rhsDissipationInviscidFine
   */
  void rhsDissipationInviscidFine();

  /**
   * \brief
   * Adds dissipation and inviscid terms to the RHS on coarse MG levels.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-18
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/rhsDissipationInviscidCoarse.C rhsDissipationInviscidCoarse
   */
  void rhsDissipationInviscidCoarse();

  /**
   * \brief
   * Adds viscous terms to the RHS on fine MG levels.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-23
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/rhsViscousFine.C rhsViscousFine
   */
  void rhsViscousFine();

  /**
   * \brief
   * Adds viscous terms to the RHS on coarse MG levels.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-23
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/rhsViscousCoarse.C rhsViscousCoarse
   */
  void rhsViscousCoarse();

  /**
   * \brief
   * Adds boundary terms to the RHS on fine MG levels.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-23
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/rhsBoundaryFine.C rhsBoundaryFine
   */
  void rhsBoundaryFine();

  /**
   * \brief
   * Adds boundary terms to the RHS on coarse MG levels.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-23
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/rhsBoundaryCoarse.C rhsBoundaryCoarse
   */
  void rhsBoundaryCoarse();

  /**
   * \brief
   * Adds physical source terms to the RHS on fine MG levels.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-23
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/rhsSourceFine.C rhsSourceFine
   */
  void rhsSourceFine();

  /**
   * \brief
   * Adds physical source terms to the RHS on coarse MG levels.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-23
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/rhsSourceCoarse.C rhsSourceCoarse
   */
  void rhsSourceCoarse();

  /**
   * \brief
   * Adds MMS source terms to the RHS on fine MG levels.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-23
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/rhsSourceMMS.C rhsSourceMMS
   */
  void rhsSourceMMS();

  /**
   * \brief
   * Adds multigrid source terms to the RHS on coarse MG levels.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-23
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Numerics/rhsSourceMG.C rhsSourceMG
   */
  void rhsSourceMG(const int& mode,
		   const int& sweep,
		   const int& linearStep);

  /**
   * \brief
   * Computes gradient limiter.
   * \param mglevel Multigrid level.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-18
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/limit.C
   */
  void limit(const int& mglevel);

  /**
   * \brief
   * Computes inviscid spectral radius.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-18
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/specRadi.C
   */
  void specRadi();

  /**
   * \brief
   * Computes viscous spectral radius.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-18
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/specRadv.C
   */
  void specRadv();

  /**
   * \brief
   * Computes source spectral radius.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-18
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/specRads.C
   */
  void specRads();

  /**
   * \brief
   * Computes LHS contributions from physical and pseudo-time terms.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-18
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/lhsTime.C
   */
  void lhsTime(const int& step,
	       const int& pseudoStep);

  /**
   * \brief
   * Computes LHS contributions from dissipation terms.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-11-2
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/lhsDissipation.C
   */
  void lhsDissipation();

  /**
   * \brief
   * Computes LHS contributions from inviscid terms.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-11-2
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/lhsInviscid.C
   */
  void lhsInviscid();

  /**
   * \brief
   * Computes LHS contributions from viscous terms on fine MG levels.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-11-2
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/lhsViscousFine.C
   */
  void lhsViscousFine();

  /**
   * \brief
   * Computes LHS contributions from viscous terms on coarse MG levels.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-11-2
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/lhsViscousCoarse.C
   */
  void lhsViscousCoarse();

  /**
   * \brief
   * Computes LHS contributions from boundary terms on coarse and
   * fine MG levels.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-11-2
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/lhsBoundary.C
   */
  void lhsBoundary(const int& mglevel);

  /**
   * \brief
   * Computes LHS contributions from source terms on fine MG levels.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-11-2
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/lhsSourceFine.C
   */
  void lhsSourceFine();

  /**
   * \brief
   * Computes LHS contributions from source terms on coarse MG levels.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-11-2
   * \par Further Documentation:
   * \par Source Code:
   * \include src/Numerics/lhsSourceCoarse.C
   */
  void lhsSourceCoarse();


  // get methods
  StrandSystem* getSystem();
  const int&    getNGfaces();
  const int&    getNFringe();
  const int&    getNPedges();
  const int&    getNEdges();
  const int&    getNqa();
  const int&    getPid();
  const int&    getNdim();
  const int&    getNBpatches();
  const int&    getInviscid();
  const int&    getViscous();
  const int&    getSource();
  const int&    getDissipation();
        int*    getFTag();
  Array1D<int>* getFTagArray();
        int*    getBTag();
  Array1D<int>* getBTagArray();
        int*    getEdge();
  Array2D<int>* getEdgeArray();
        double* getFacu();
  Array3D<double>* getFacuArray();
        double* getFacs();
  Array3D<double>* getFacsArray();
        double* getXvu();
  Array2D<double>* getXvuArray();
        double* getXvs();
  Array2D<double>* getXvsArray();
        double* getR();
  Array3D<double>* getRArray();
        double* getQa();
  Array3D<double>* getQaArray();
  int*          getLimFlag();
  int*          getNodalQFlag();
  int*          getNodalQaFlag();
  int*          getGradQFlag();
  int*          getGradQaFlag();


  // index functions to deal with multidimensional arrays
  void indlsp(const int& k, const int& j, const int& n, int& i);


  // strand solver numerics inputs
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
  static int    nLevels;
  static int    mgCycle;
  static int    gradient;
  static int    nodeVal;
  static int    perturb;
  static int    limiter;
  static int    nRamp;
  static double brelax;
  static double gradClip;
  static double dtUnsteady;
  static double cflLinear;
  static double cfl0;
  static double vnn0;
  static double cfl;
  static double vnn;
  static double convLimit;
  static double relax;
  static double coarseDis;

  // dimensions for local solver data
  int     nq;
  int     nqa;
  int     ndim;
  int     nBpatches;
  int     nDofs;
  int     nFaces;
  int     nGfaces;
  int     nNodes;
  int     nGnodes;
  int     nEdges;
  int     nBedges;
  int     nPedges;
  int     nPstr;
  int     nFringe;
  int     nSharp;
  int     nFacesP;
  int     nGfacesP;
  int     nEdgesP;
  int     nBedgesP;
  int     nPedgesP;
  int     nPstrP;
  int     nFringeP;

  // local solver data
  StrandSystem* sys;
  int     pid;
  int     inviscid;
  int     dissipation;
  int     viscous;
  int     source;
  int     sourceMMS;
  int     gradQFlag;
  int     gradQaFlag;
  int     nodalQFlag;
  int     nodalQaFlag;
  int     limFlag;
  int*    gradQFlagP;
  int*    gradQaFlagP;
  int*    nodalQFlagP;
  int*    nodalQaFlagP;
  int*    limFlagP;
  int     nqGradQ;
  int     nqaGradQa;
  double  forcex;
  double  forcey;
  Array1D<double> dlim;
  Array1D<int> iqgrad;
  Array1D<int> iqagrad;
  Array1D<double> rmsNorm;
  Array1D<double> rms;

  Array2D<int> face;
  Array1D<int> fTag;
  Array1D<int> bTag;
  Array1D<int> fClip;
  Array1D<int> sFlag;
  Array2D<int> edge;
  Array2D<int> edgp;
  Array1D<int> edgn;
  Array1D<int> ncsc;
  Array1D<int> csc;
  Array1D<int> ncsp;
  int**   csp;
  Array1D<int> gsMap;

  Array1D<int> f2cc;
  Array1D<int> f2ce;
  Array1D<int> f2cs;

  Array3D<double> x;
  Array3D<double> x0;
  Array3D<double> x1;
  Array3D<double> x2;
  Array3D<double> xc;
  Array2D<double> v;
  Array2D<double> v1;
  Array2D<double> v2;
  Array3D<double> facu;
  Array3D<double> facs;
  double** lsp;
  Array1D<double> xStr;

  Array2D<double> xvu;
  Array2D<double> xvs;
  Array2D<double> dvu;
  Array2D<double> dvs;
  Array3D<double> nvu;
  Array3D<double> nvs;

  Array3D<double> q;
  Array3D<double> qa;
  Array3D<double> q0;
  Array3D<double> q1;
  Array3D<double> q2;
  Array3D<double> qp;
  Array3D<double> qap;
  Array4D<double> qx;
  Array4D<double> qax;
  Array3D<double> r;
  Array3D<double> dq;
  Array3D<double> s;
  Array2D<double> dt;
  Array2D<double> radi;
  Array2D<double> radv;
  Array2D<double> rads;
  Array4D<double> dd;
  Array4D<double> dm;
  Array4D<double> dp;
  Array4D<double> bu;
  Array3D<double> fwc;
  Array3D<double> limu;
  Array3D<double> lims;

  const Array2D<int>* edgeP;
  const Array1D<int>* bTagP;
  const Array1D<int>* fTagP;
  const Array1D<int>* fClipP;
  const Array3D<double>* facuP;
  const Array3D<double>* facsP;
  const Array2D<double>* xvuP;
  const Array2D<double>* xvsP;
        Array3D<double>* qP;
        Array3D<double>* qaP;
  const Array3D<double>* rP;
  const Array2D<double>* vP;
};
#endif
