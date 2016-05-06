/**
 * \brief
 * Class Strand2dFCSystem holds the data and specifies the operations for
 * various systems of equations, such as single phase, turbulent, or
 * multi-species reacting gas.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-11
 */


#ifndef included_Strand2dFCSystem
#define included_Strand2dFCSystem

#include "STRAND2DFC_defs.h"


class Strand2dFCSystem
{
 public:

  /**
   * \brief
   * Constructor for the Strand2dFCSystem class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Strand2dFCSystem/Strand2dFCSystem.C Strand2dFCSystem
   */
  Strand2dFCSystem();

  /**
   * \brief
   * Destructor for the Strand2dFCSystem class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Strand2dFCSystem/Strand2dFCSystem.C ~Strand2dFCSystem
   */
  virtual ~Strand2dFCSystem();

  /**
   * \brief
   * Reads inputs for the Strand2dFCSystem layer.
   * \param inputFile Name of input file.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-01-02
   * \par Further Documentation:
   */
  virtual void inputRead(const string& inputFile)=0;

  /**
   * \brief
   * Provides relevant StrandSystem data to the Numerics routines
   * \param iPrint Specifes whether inputs are printed out; usually iprint
   * is set to 1 in node 0.
   * \param iTest Specifies whether run is in unit test mode.
   * \param iDebug Specifies whether run is in debug mode, i.e., triggers
   * more verbose data.
   * \param tmp Temporary dimension of iqgrad,iqgrada. Typically set to
   * something large, like 500. This is needed because the main program
   * does not know the dimension of iqgrad,iqagrad yet.
   * \param nq Total number of equations or Q variables in system.
   * \param nqa Total number of additional variables to be stored.
   * \param ndim Number of spatial dimensions.
   * \param inviscid Flag to add inviscid terms (0 do not add, 1 add).
   * \param viscous Flag to add viscous terms (0 do not add, 1 add).
   * \param source Flag to add physical source terms
   * (0 do not add, 1 physical source).
   * \param sourceMMS Flag to add MMS source terms
   * (0 do not add, 1 MMS source).
   * \param dissipation Flag to add dissipation (0 do not add, 1 add).
   * \param nBPatches Number of boundary patches.
   * \param iqgrad Dimension (tmp), flag for whether gradient of Q is
   * required: 0 or 1 for each.
   * \param iqagrad Dimension (tmp), flag for whether gradient of Qa is
   * required: 0 or 1 for each.
   * \param dlim Constant used in limiter computations.
   * \param dqNorm Normalization values for Q.
   * \param nOutputVars Number of output variables.
   * \param outputVars List of output variables.
   * \param outputVarLength List of output variable lengths (scalar or vector).
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/StrandSystem/StrandSystem.C prepSetup
   */
  virtual void prepSetup(const int& iPrint,
                         const int& iTest,
                         const int& iDebug,
                         const int& tmp,
                         const int& iSolnFile,
                         const int& iResdFile,
                         const int& iErrFile,
                         int& nq,
                         int& nqa,
                         int& ndim,
                         int& inviscid,
                         int& viscous,
                         int& source,
                         int& sourceMMS,
                         int& dissipation,
                         int& nBpatches,
                         int* iqgrad,
                         int* iqagrad,
                         double* dlim,
                         double* rmsNorm,
			 int& nOutputVars,
			 int* outputVarLength,
			 string* outputVars)=0;

  /**
   * \brief
   * Initializes q over a set of dof locations.
   * \param npts Number of points at which to set q.
   * \param x Coordinates of each dof, laid out in (x,y) pairs.
   * \param q Solution vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-25
   * \par Further Documentation:
   */
  virtual void initFlow(const int& npts,
                        const double* x,
                        double* q)=0;

  /**
   * \brief
   * Initializes source terms over a set of dof locations.
   * \param npts Number of points at which to set q.
   * \param x Coordinates of each dof, laid out in (x,y) pairs.
   * \param s Source term vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-25
   * \par Further Documentation:
   */
  virtual void initSource(const int& npts,
                          const double* x,
                          double* s)=0;

  /**
   * \brief
   * Initializes viscous fluxes.
   * \param npts Number of points at which to set q.
   * \param x Coordinates of each dof, laid out in (x,y) pairs.
   * \param f Viscous flux vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-25
   * \par Further Documentation:
   */
  virtual void initVisFlux(const int& npts,
			   const double* x,
			   double* f)=0;

  /**
   * \brief
   * Initializes viscous fluxes.
   * \param npts Number of points at which to set q.
   * \param tag Boundary tag.
   * \param x Coordinates of each dof, laid out in (x,y) pairs.
   * \param f Solution vector at boundary.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-25
   * \par Further Documentation:
   */
  virtual void initPenaltyData(const int& npts,
			       const int* tag,
			       const double* x,
			       double* f)=0;

  /**
   * \brief
   * Computes additional variables as a function of Q.
   * \param npts Number of points at which to set qa.
   * \param q Q vector.
   * \param qa Additional Q vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-10
   * \par Further Documentation:
   */
  virtual void stepQAdd(const int& npts,
                        const double* q,
                        double* qa)=0;

  /**
   * \brief
   * Computes the X-inviscid flux.
   * \param npts Number of points at which to compute the flux.
   * \param q Q vector.
   * \param qa Additional Q vector.
   * \param f X-inviscid flux vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-10
   * \par Further Documentation:
   */
  virtual void rhsInvFluxX(const int& npts,
			   const double* q,
			   const double* qa,
			   double* f)=0;

  /**
   * \brief
   * Computes the Y-inviscid flux.
   * \param npts Number of points at which to compute the flux.
   * \param q Q vector.
   * \param qa Additional Q vector.
   * \param g Y-inviscid flux vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-10
   * \par Further Documentation:
   */
  virtual void rhsInvFluxY(const int& npts,
			   const double* q,
			   const double* qa,
			   double* g)=0;

  /**
   * \brief
   * Computes the inviscid flux Jacobian.
   * \param npts Number of points at which to compute the Jacobian.
   * \param nx X-component of direction vector.
   * \param ny Y-component of direction vector.
   * \param q Q vector.
   * \param qa Additional Q vector.
   * \param A inviscid flux Jacobian.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2014-6-16
   * \par Further Documentation:
   */
  virtual void lhsInvFluxJacobian(const int& npts,
				  const double* nx,
				  const double* ny,
				  const double* q,
				  const double* qa,
				  double* A)=0;

  /**
   * \brief
   * Computes the Jacobian of the viscous gradient variables with
   * respect to the conserved variables.
   * \param npts Number of points at which to compute the Jacobian.
   * \param q Q vector.
   * \param qa Additional Q vector.
   * \param A Viscous variable Jacobian.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2014-6-23
   * \par Further Documentation:
   */
  virtual void lhsVisVariableJacobian(const int& npts,
				      const double* q,
				      const double* qa,
				      double* A)=0;

  /**
   * \brief
   * Computes the X/Y-viscous flux.
   * \param npts Number of points at which to compute the flux.
   * \param q Q vector.
   * \param qa Additional Q vector.
   * \param qax dQa/dx vector.
   * \param qay dQa/dy vector.
   * \param f X-viscous flux vector.
   * \param g Y-viscous flux vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-10
   * \par Further Documentation:
   */
  virtual void rhsVisFlux(const int& npts,
			  const double* q,
			  const double* qa,
			  const double* qax,
			  const double* qay,
			  double* f,
			  double* g)=0;

  /**
   * \brief
   * Computes the X-viscous flux.
   * \param npts Number of points at which to compute the flux.
   * \param q Q vector.
   * \param qa Additional Q vector.
   * \param qax dQa/dx vector.
   * \param qay dQa/dy vector.
   * \param f X-viscous flux vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-10
   * \par Further Documentation:
   */
  virtual void rhsVisFluxX(const int& npts,
			   const double* q,
			   const double* qa,
			   const double* qax,
			   const double* qay,
			   double* f)=0;

  /**
   * \brief
   * Computes the Y-viscous flux.
   * \param npts Number of points at which to compute the flux.
   * \param q Q vector.
   * \param qa Additional Q vector.
   * \param qax dQa/dx vector.
   * \param qay dQa/dy vector.
   * \param g Y-viscous flux vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-10
   * \par Further Documentation:
   */
  virtual void rhsVisFluxY(const int& npts,
			   const double* q,
			   const double* qa,
			   const double* qax,
			   const double* qay,
			   double* g)=0;

  /**
   * \brief
   * Computes the s-direction viscous flux component.
   * \param npts Number of points at which to compute the flux.
   * \param jac Local jacobian.
   * \param xs Local dx/ds.
   * \param ys Local dy/ds.
   * \param xn Local dx/dn.
   * \param yn Local dy/dn.
   * \param q Q vector.
   * \param qa Additional Q vector.
   * \param qas dQa/ds vector.
   * \param f s-viscous flux vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-12-3
   * \par Further Documentation:
   */
  virtual void rhsVisFluxS(const int& npts,
			   const double* jac,
			   const double* xs,
			   const double* ys,
			   const double* xn,
			   const double* yn,
			   const double* q,
			   const double* qa,
			   const double* qas,
			   double* f)=0;

  /**
   * \brief
   * Computes the n-direction viscous flux coefficient matrix Jacobian.
   * \param npts Number of points at which to compute the matrix Jacobian.
   * \param jac Local jacobian.
   * \param xs Local dx/ds.
   * \param ys Local dy/ds.
   * \param xn Local dx/dn.
   * \param yn Local dy/dn.
   * \param q Q vector at Jacobian evaluation point.
   * \param qa Additional Q vector at Jacobian evaluation point.
   * \param qi Q vector at product evaluation point.
   * \param qai Additional Q vector at product evaluation point.
   * \param A n-viscous flux coefficient matrix Jacobian.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2014-6-23
   * \par Further Documentation:
   */
  virtual void lhsVisFluxNCoeffJacobian(const int& npts,
					const double* jac,
					const double* xs,
					const double* ys,
					const double* xn,
					const double* yn,
					const double* q,
					const double* qa,
					const double* qi,
					const double* qai,
					double* A)=0;

  /**
   * \brief
   * Computes the n-direction viscous flux coefficient matrix.
   * \param npts Number of points at which to compute the matrix.
   * \param jac Local jacobian.
   * \param xs Local dx/ds.
   * \param ys Local dy/ds.
   * \param xn Local dx/dn.
   * \param yn Local dy/dn.
   * \param q Q vector.
   * \param qa Additional Q vector.
   * \param b n-viscous flux coefficient matrix.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-12-3
   * \par Further Documentation:
   */
  virtual void rhsVisFluxNCoeff(const int& npts,
				const double* jac,
				const double* xs,
				const double* ys,
				const double* xn,
				const double* yn,
				const double* q,
				const double* qa,
				double* b)=0;

  /**
   * \brief
   * Computes directed dissipation flux vector.
   * \param npts Number of points at which to compute the flux.
   * \param nx X-component of direction vector.
   * \param ny Y-component of direction vector.
   * \param qL Left Q vector.
   * \param qR Right Q vector.
   * \param dq qR-qL for a given edge.
   * \param f directed dissipation flux vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-10
   * \par Further Documentation:
   */
  virtual void rhsDisFlux(const int& npts,
			  const double* nx,
			  const double* ny,
			  const double* qL,
			  const double* qR,
			  const double* dq,
			  double* f)=0;

  /**
   * \brief
   * Computes physical source terms.
   * \param npts Number of points at which to compute the source.
   * \param q State vector.
   * \param qa Additional state vector.
   * \param qax X-component of gradient of qa.
   * \param qay Y-component of gradient of qa.
   * \param s source vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2014-6-26
   * \par Further Documentation:
   */
  virtual void rhsSource(const int& npts,
			 const double* q,
			 const double* qa,
			 const double* qax,
			 const double* qay,
			 double* s);

  /**
   * \brief
   * Computes physical source term Jacobians.
   * \param npts Number of points at which to compute the source Jacobian.
   * \param q State vector.
   * \param qa Additional state vector.
   * \param qax X-component of gradient of qa.
   * \param qay Y-component of gradient of qa.
   * \param A source vector Jacobian.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2014-6-26
   * \par Further Documentation:
   */
  virtual void lhsSourceJacobian(const int& npts,
				 const double* q,
				 const double* qa,
				 const double* qax,
				 const double* qay,
				 double* A);

  /**
   * \brief
   * Sets wall distance and stores in qa.
   * \param npts Number of points at which to set qa.
   * \param dw Distance to the nearest wall.
   * \param qa Additional Q vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2014-6-24
   * \par Further Documentation:
   */
  virtual void initWallDist(const int& npts,
                            const double* dw,
                            double* qa);


  /**
   * \brief
   * Computes directed dissipation flux vector Jacobian.
   * \param npts Number of points at which to compute the flux.
   * \param nx X-component of direction vector.
   * \param ny Y-component of direction vector.
   * \param qL Left Q vector.
   * \param qR Right Q vector.
   * \param A directed dissipation flux vector Jacobian.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2014-06-16
   * \par Further Documentation:
   */
  virtual void lhsDisFluxJacobian(const int& npts,
				  const double* nx,
				  const double* ny,
				  const double* qL,
				  const double* qR,
				  double* A)=0;

  /**
   * \brief
   * Computes invsicid spectral radius.
   * \param npts Number of points at which to compute the flux.
   * \param aL Weight coefficient for left state.
   * \param aR Weight coefficient for right state.
   * \param Ax X-component of area vector.
   * \param Ay Y-component of area vector.
   * \param qL Left Q vector.
   * \param qR Right Q vector.
   * \param qaL Left Qa vector.
   * \param qaR Right Qa vector.
   * \param sr Inviscid spectral radius.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-10
   * \par Further Documentation:
   */
  virtual void stepInvEigenvalue(const int& npts,
				 const double& aL,
				 const double& aR,
				 const double* Ax,
				 const double* Ay,
				 const double* qL,
				 const double* qR,
				 const double* qaL,
				 const double* qaR,
				 double* sr)=0;

  /**
   * \brief
   * Computes viscous spectral radius.
   * \param npts Number of points at which to compute the flux.
   * \param aL Weight coefficient for left state.
   * \param aR Weight coefficient for right state.
   * \param Ax X-component of area vector.
   * \param Ay Y-component of area vector.
   * \param qL Left Q vector.
   * \param qR Right Q vector.
   * \param qaL Left Qa vector.
   * \param qaR Right Qa vector.
   * \param vL Left volume.
   * \param vR Right volume.
   * \param sr Inviscid spectral radius.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-10
   * \par Further Documentation:
   */
  virtual void stepVisEigenvalue(const int& npts,
				 const double& aL,
				 const double& aR,
				 const double* Ax,
				 const double* Ay,
				 const double* qL,
				 const double* qR,
				 const double* qaL,
				 const double* qaR,
				 const double* vL,
				 const double* vR,
				 double* sr)=0;

  /**
   * \brief
   * Resolves boundary conflicts at sharp corners.
   * \param npts Number of boundary points to update.
   * \param tag Current boundary tag.
   * \param tagC Candidate boundary tag.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-04-24
   * \par Further Documentation:
   */
  virtual void bcConflict(const int& npts,
			  int* tag,
			  const int* tagC)=0;

  /**
   * \brief
   * Provides the BC Vector (boundary residual).
   * (Neumann gradient conditions to be added soon)
   * \param npts Number of boundary points to update.
   * \param tag Boundary tag.
   * \param nx Outward boundary normal in (nx,ny) pairs.
   * \param q Q vector at boundary.
   * \param qa Qa vector at boundary.
   * \param rb boundary residual vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-04-24
   * \par Further Documentation:
   */
  virtual void rhsBCVector(const int& npts,
			   const int* tag,
			   const double* nx,
			   const double* q,
			   const double* qa,
			   double* rb)=0;

  /**
   * \brief
   * Provides the Jacobian of the BC Vector with respect to Q at the boundary.
   * (Neumann gradient conditions to be added soon)
   * \param npts Number of boundary points to update.
   * \param tag Boundary tag.
   * \param nx Outward boundary normal in (nx,ny) pairs.
   * \param q Q vector at boundary.
   * \param qa Qa vector at boundary.
   * \param M boundary residual vector Jacobian.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-04-24
   * \par Further Documentation:
   */
  virtual void lhsBCVectorSelfJacobian(const int& npts,
				       const int* tag,
				       const double* nx,
				       const double* q,
				       const double* qa,
				       double* M)=0;

  /**
   * \brief
   * Provides the BC selection matrix.
   * \param npts Number of boundary points to update.
   * \param tag Boundary tag.
   * \param nx Outward boundary normal in (nx,ny) pairs.
   * \param q Q vector at boundary.
   * \param qa Qa vector at boundary.
   * \param L Selection matrix.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-04-24
   * \par Further Documentation:
   */
  virtual void rhsBCSelectionMatrix(const int& npts,
				    const int* tag,
				    const double* nx,
				    const double* q,
				    const double* qa,
				    double* L)=0;

  /**
   * \brief
   * Provides the BC penalty term.
   * \param npts Number of boundary points to update.
   * \param inout Flag to determine incoming or outgoing characteristics.
   * \param tag Boundary tag.
   * \param A Directed flux area components in (Ax,Ay) pairs.
   * \param Pinv0 Inverse of norm at boundary.
   * \param q Q vector at boundary.
   * \param qa Qa vector at boundary.
   * \param g Boundary data.
   * \param uw Wall velocity.
   * \param rb boundary penalty term.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2014-6-2
   * \par Further Documentation:
   */
  virtual void rhsBCPenalty(const int& npts,
			    const int* tag,
			    const int& inout,
			    const double* A,
			    const double* Pinv0,
			    const double* q,
			    const double* qa,
			    const double* g,
			    const double* uw,
			    double* rb)=0;

  /**
   * \brief
   * Provides the BC penalty term.
   * \param npts Number of boundary points to update.
   * \param tag Boundary tag.
   * \param Pinv0 Inverse of norm at boundary.
   * \param q Q vector at boundary.
   * \param qa Qa vector at boundary.
   * \param gv Viscous flux penalty.
   * \param uw Wall velocity.
   * \param rb boundary penalty term.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2014-6-2
   * \par Further Documentation:
   */
  virtual void rhsBCPenaltyVis(const int& npts,
			       const int* tag,
			       const double* Pinv0,
			       const double* q,
			       const double* qa,
			       const double* gv,
			       const double* uw,
			       double* rb)=0;

  /**
   * \brief
   * Provides the BC penalty term Jacobian.
   * \param npts Number of boundary points to update.
   * \param inout Flag to determine incoming or outgoing characteristics.
   * \param tag Boundary tag.
   * \param A Directed flux area components in (Ax,Ay) pairs.
   * \param Pinv0 Inverse of norm at boundary.
   * \param q Q vector at boundary.
   * \param qa Qa vector at boundary.
   * \param g Boundary data.
   * \param uw Wall velocity.
   * \param M boundary penalty term Jacobian.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2014-6-2
   * \par Further Documentation:
   */
  virtual void lhsBCPenaltyJacobian(const int& npts,
				    const int* tag,
				    const int& inout,
				    const double* A,
				    const double* Pinv0,
				    const double* q,
				    const double* qa,
				    const double* g,
				    const double* uw,
				    double* M)=0;

  /**
   * \brief
   * Prints out a block of specified output data.
   * npts Number of output points.
   * ffile Ouput file stream.
   * outputVar String containing variable to output.
   * q Q vector.
   * qa Qa vector.
   * qe Exact Q vector.
   * r Residual vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-10
   * \par Further Documentation:
   */
  virtual void outputSolution(const int& npts,
			      ofstream& ffile,
			      const string& outputVar,
			      const double* q,
			      const double* qa,
			      const double* e,
			      const double* r)=0;

  /**
   * \brief
   * Returns surface forces.
   * \param npts Number of surface points at which to compute forces.
   * \param tag Boundary tag.
   * \param xs dx/ds surface mapping.
   * \param ys dy/ds surface mapping.
   * \param q Q vector.
   * \param qa Qa vector.
   * \param qx X-component of gradient of Q vector.
   * \param qy Y-component of gradient of Q vector.
   * \param qax X-component of gradient of Qa vector.
   * \param qay Y-component of gradient of Qa vector.
   * \param force Force vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-5-28
   * \par Further Documentation:
   */
  virtual void outputSurfaceForces(const int& npts,
				   const int* tag,
				   const double* xs,
				   const double* ys,
				   const double* q,
				   const double* qa,
				   const double* qx,
				   const double* qy,
				   const double* qax,
				   const double* qay,
				   double* force)=0;

  /**
   * \brief
   * Finalizes the Strand2dFCSystem class instance.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-1-1
   * \par Further Documentation:
   */
  virtual void finalize()=0;


 protected:

  int iPrint;
  int iTest;
  int iDebug;
  int iSolnFile;
  int iResdFile;
  int iErrFile;
  int nq;
  int nqa;
  int nqGradQ;
  int nqaGradQa;
  int ndim;
  int ncomp;
  int isolution;
  int nBpatches;
  int inviscid;
  int viscous;
  int source;
  int sourceMMS;
  int dissipation;


 private:


};
#endif
