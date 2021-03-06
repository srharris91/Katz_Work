/**
 * \brief
 * Class StrandSystem holds the data and specifies the operations for various
 * systems of equations, such as single phase, turbulent, or multi-species
 * reacting gas.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-11
 */


#ifndef included_StrandSystem
#define included_StrandSystem

#include "STRAND_defs.h"


class StrandSystem
{
 public:

  /**
   * \brief
   * Constructor for the StrandSystem class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/StrandSystem/StrandSystem.C StrandSystem
   */
  StrandSystem();

  /**
   * \brief
   * Destructor for the StrandSystem class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/StrandSystem/StrandSystem.C ~StrandSystem
   */
  virtual ~StrandSystem();

  /**
   * \brief
   * Reads inputs for the StrandSystem layer.
   * \param inputFile name of StrandSystem input file.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/StrandSystem/StrandSystem.C inputRead
   */
  virtual void inputRead(const string& inputFile);

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
   * \param dqNorm normalization values for Q.
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
			 int& nq,
			 int& nqa,
			 int& ndim,
			 int& inviscid,
			 int& viscous,
			 int& source,
			 int& sourceMMS,
			 int& dissipation,
			 int& nBpatches,
			 int* iqgradT,
			 int* iqagradT,
			 double* dlim,
			 double* rmsNorm);

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
   * Initializes source terms over a set of boundary dof locations.
   * \param npts Number of points at which to set q.
   * \param tag Tag number of the boundary dofs.
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
			  const int* tag,
			  const double* x,
			  double* s)=0;

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
   * 2012-10-10
   * \par Further Documentation:
   */
  virtual void initWallDist(const int& npts,
			    const double* dw,
			    double* qa);

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
   * Provides max inviscid eigenvalue at an edge separating volumes 1 and 2.
   * \param npts Number of edges at which to compute spectral radius.
   * \param A Face area in (Ax,Ay) pairs.
   * \param xv Mesh velocity dotted with face area.
   * \param q1 Q vector at volume 1.
   * \param q2 Q vector at volume 2.
   * \param qa1 Qa vector at volume 1.
   * \param qa2 Qa vector at volume 2.
   * \param sr Spectral radius.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-10
   * \par Further Documentation:
   */
  virtual void stepInvEigenvalue(const int& npts,
				 const double* A,
				 const double* xv,
				 const double* q1,
				 const double* q2,
				 const double* qa1,
				 const double* qa2,
				 double* sr)=0;

  /**
   * \brief
   * Provides max viscous eigenvalue at an edge separating volumes 1 and 2.
   * \param npts Number of edges at which to compute spectral radius.
   * \param A Face area in (Ax,Ay) pairs.
   * \param q1 Q vector at volume 1.
   * \param q2 Q vector at volume 2.
   * \param qa1 Qa vector at volume 1.
   * \param qa2 Qa vector at volume 2.
   * \param v1 Volume 1.
   * \param v2 Volume 2.
   * \param sr Spectral radius.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-10
   * \par Further Documentation:
   */
  virtual void stepVisEigenvalue(const int& npts,
				 const double* A,
				 const double* q1,
				 const double* q2,
				 const double* qa1,
				 const double* qa2,
				 const double* v1,
				 const double* v2,
				 double* sr)=0;

  /**
   * \brief
   * Provides max source eigenvalue.
   * \param npts Number of points at which to compute spectral radius.
   * \param q Q vector.
   * \param qa Qa vector.
   * \param sr Spectral radius.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-10
   * \par Further Documentation:
   */
  virtual void stepSourceEigenvalue(const int& npts,
				    const double* q,
				    const double* qa,
				    double* sr)=0;

  /**
   * \brief
   * Provides inviscid flux at interface given a series of left and right states.
   * \param npts Number of points at which to compute spectral radius.
   * \param A Face area in (Ax,Ay) pairs.
   * \param xv Mesh velocity dotted with face area.
   * \param ql Q vector at left state.
   * \param qr Q vector at right state.
   * \param f Inviscid flux vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-10
   * \par Further Documentation:
   */
  virtual void rhsInvFlux(const int& npts,
			  const double* A,
			  const double* xv,
			  const double* ql,
			  const double* qr,
			  double* f)=0;

  /**
   * \brief
   * Provides disspation flux at interface given a series of left and
   * right states.
   * \param npts Number of points at which to compute spectral radius.
   * \param A Face area in (Ax,Ay) pairs.
   * \param xv Mesh velocity dotted with face area.
   * \param ql Q vector at left state.
   * \param qr Q vector at right state.
   * \param f Dissipation flux vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-10
   * \par Further Documentation:
   */
  virtual void rhsDisFlux(const int& npts,
			  const double* A,
			  const double* xv,
			  const double* ql,
			  const double* qr,
			  double* f)=0;

  /**
   * \brief
   * Provides disspation flux at interface given a series of left and
   * right states on coarse MG levels.
   * \param npts Number of points at which to compute spectral radius.
   * \param A Face area in (Ax,Ay) pairs.
   * \param xv Mesh velocity dotted with face area.
   * \param ql Q vector at left state.
   * \param qr Q vector at right state.
   * \param f Dissipation flux vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-10
   * \par Further Documentation:
   */
  virtual void rhsDisFluxCoarse(const int& npts,
				const double* A,
				const double* xv,
				const double* ql,
				const double* qr,
				double* f)=0;

  /**
   * \brief
   * Provides viscous flux at interface given interface states and gradients.
   * \param npts Number of points at which to compute spectral radius.
   * \param A Face area in (Ax,Ay) pairs.
   * \param q Q vector at interface.
   * \param qa Qa vector at interface.
   * \param qx X-component of gradient of Q vector at interface.
   * \param qy Y-component of gradient of Q vector at interface.
   * \param qax X-component of gradient of Qa vector at interface.
   * \param qay Y-component of gradient of Qa vector at interface.
   * \param f Viscous flux vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-23
   * \par Further Documentation:
   */
  virtual void rhsVisFlux(const int& npts,
			  const double* A,
			  const double* q,
			  const double* qa,
			  const double* qx,
			  const double* qy,
			  const double* qax,
			  const double* qay,
			  double* f)=0;

  /**
   * \brief
   * Provides viscous flux at interface on coarse MG levels.
   * \param npts Number of points at which to compute spectral radius.
   * \param A Face area in (Ax,Ay) pairs.
   * \param vl Volume at left state.
   * \param vr Volume at right state.
   * \param ql Q vector at left state.
   * \param qr Q vector at right state.
   * \param qal Qa vector at left state.
   * \param qar Qa vector at right state.
   * \param f Viscous flux vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-23
   * \par Further Documentation:
   */
  virtual void rhsVisFluxCoarse(const int& npts,
				const double* A,
				const double* vl,
				const double* vr,
				const double* ql,
				const double* qr,
				const double* qal,
				const double* qar,
				double* f)=0;

  /**
   * \brief
   * Provides BC vector for a series of Q and extrapolated values 
   * (Neumann gradient conditions to be added soon)
   * \param npts Number of points at which to compute spectral radius.
   * \param nx Outward boundary normal in (nx,ny) pairs.
   * \param wx Wall velocity in (wx,wy) pairs.
   * \param qe Extrapolated Q vector at boundary location.
   * \param qae Extrapolated Qa vector at boundary location.
   * \param q Q vector of boundary dof.
   * \param qa Qa vector of boundary dof.
   * \param r BC vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-27
   * \par Further Documentation:
   */
  virtual void rhsBCVector(const int& npts,
			   const int* tag,
			   const double* nx,
			   const double* wx,
			   const double* qe,
			   const double* qae,
			   const double* q,
			   const double* qa,
			   double* r)=0;

  /**
   * \brief
   * Provides RHS source term computation
   * \param npts Number of points at which to compute spectral radius.
   * \param q Q vector.
   * \param qa Qa vector.
   * \param qx Qx vector containing x and y gradient components.
   * \param qax Qax vector containing x and y gradient components.
   * \param f source term vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-27
   * \par Further Documentation:
   */
  virtual void rhsSource(const int& npts,
			 const double* q,
			 const double* qa,
			 const double* qx,
			 const double* qax,
			 double* r);

  /**
   * \brief
   * Provides RHS source term computation on coarse MG levels.
   * \param npts Number of points at which to compute spectral radius.
   * \param v Cell volume.
   * \param q Q vector.
   * \param qa Qa vector.
   * \param f source term vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-27
   * \par Further Documentation:
   */
  virtual void rhsSourceCoarse(const int& npts,
			       const double* v,
			       const double* q,
			       const double* qa,
			       double* r);

  /**
   * \brief
   * Provides inviscid flux Jacobian for a series of locations.
   * \param npts Number of points at which to compute inviscid Jacobian.
   * \param A Face area in (Ax,Ay) pairs.
   * \param xv Mesh velocity dotted with face area.
   * \param q Q vector.
   * \param qa Qa vector.
   * \param M Inviscid flux Jacobian matrix (row order).
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-30
   * \par Further Documentation:
   */
  virtual void lhsInvFluxJacobian(const int& npts,
				  const double* A,
				  const double* xv,
				  const double* q,
				  const double* qa,
				  double* M)=0;

  /**
   * \brief
   * Provides dissipative flux Jacobian for a series of locations.
   * \param npts Number of points at which to compute dissipative Jacobian.
   * \param A Face area in (Ax,Ay) pairs.
   * \param xv Mesh velocity dotted with face area.
   * \param ql Q vector at left state.
   * \param qr Q vector at right state.
   * \param M Inviscid flux Jacobian matrix (row order).
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-30
   * \par Further Documentation:
   */
  virtual void lhsDisFluxJacobian(const int& npts,
				  const double* A,
				  const double* xv,
				  const double* ql,
				  const double* qr,
				  double* M)=0;

  /**
   * \brief
   * Provides viscous flux Jacobian at interface.
   * \param npts Number of points at which to compute viscous Jacobian.
   * \param A Face area in (Ax,Ay) pairs.
   * \param B Vector with which to dot x,y viscous matrices in (Bx,By) pairs.
   * \param qe Q vector at cell interface.
   * \param qae Qa vector at cell interface.
   * \param M Viscous flux Jacobian matrix (row order).
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-27
   * \par Further Documentation:
   */
  virtual void lhsVisFluxJacobian(const int& npts,
				  const double* A,
				  const double* B,
				  const double* qe,
				  const double* qae,
				  double* M)=0;

  /**
   * \brief
   * Provides viscous flux Jacobian at interface.
   * \param npts Number of points at which to compute spectral radius.
   * \param A Face area in (Ax,Ay) pairs.
   * \param q Q vector at interface.
   * \param qa Qa vector at interface.
   * \param qx X-component of gradient of Q vector at interface.
   * \param qy Y-component of gradient of Q vector at interface.
   * \param qax X-component of gradient of Qa vector at interface.
   * \param qay Y-component of gradient of Qa vector at interface.
   * \param f Viscous flux vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-27
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/StrandSystem/StrandSystem.C lhsVisFluxJacobian
   */
  virtual void lhsVisFluxJacobian(const int& npts,
				  const double* A,
				  const double* q,
				  const double* qa,
				  const double* qx,
				  const double* qy,
				  const double* qax,
				  const double* qay,
				  double* M);

  /**
   * \brief
   * Provides source term Jacobian at interface.
   * \param npts Number of points at which to compute source Jacobian.
   * \param v Cell Volume.
   * \param q Q vector at cell interface.
   * \param qa Qa vector at cell interface.
   * \param qx Qx vector containing x and y gradient components.
   * \param qax Qax vector containing x and y gradient components.
   * \param M Viscous flux Jacobian matrix (row order).
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-27
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/StrandSystem/StrandSystem.C lhsSourceJacobian
   */
  virtual void lhsSourceJacobian(const int& npts,
				 const double* v,
				 const double* q,
				 const double* qa,
				 const double* qx,
				 const double* qax,
				 double* M);

  /**
   * \brief
   * Provides LHS source Jacobian computation on coarse MG levels.
   * \param npts Number of points at which to compute spectral radius.
   * \param v Cell volume.
   * \param q Q vector.
   * \param qa Qa vector.
   * \param f source term vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-27
   * \par Further Documentation:
   */
  virtual void lhsSourceJacobianCoarse(const int& npts,
				       const double* v,
				       const double* q,
				       const double* qa,
				       double* M);

  /**
   * \brief
   * Provides boundary Jacobian contribution for the boundary dof itself.
   * \param npts Number of points at which to compute boundary Jacobian.
   * \param nx Outward boundary normal in (nx,ny) pairs.
   * \param qe Extrapolated Q vector at boundary location.
   * \param qae Extrapolated Qa vector at boundary location.
   * \param q Q vector of boundary dof.
   * \param qa Qa vector of boundary dof.
   * \param M Boundary condition Jacobian matrix (row order) for self.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-27
   * \par Further Documentation:
   */
  virtual void lhsBCVectorSelfJacobian(const int& npts,
				       const int* tag,
				       const double* nx,
				       const double* qe,
				       const double* qae,
				       const double* q,
				       const double* qa,
				       double* M)=0;

  /**
   * \brief
   * Provides boundary Jacobian contribution for the interior location.
   * \param npts Number of points at which to compute boundary Jacobian.
   * \param nx Outward boundary normal in (nx,ny) pairs.
   * \param qe Extrapolated Q vector at boundary location.
   * \param qae Extrapolated Qa vector at boundary location.
   * \param q Q vector of boundary dof.
   * \param qa Qa vector of boundary dof.
   * \param M Boundary condition Jacobian matrix (row order) for interior.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-27
   * \par Further Documentation:
   */
  virtual void lhsBCVectorInteriorJacobian(const int& npts,
					   const int* tag,
					   const double* nx,
					   const double* qe,
					   const double* qae,
					   const double* q,
					   const double* qa,
					   double* M)=0;

  /**
   * \brief
   * Provides preconditioning Jacobian.
   * \param npts Number of points at which to compute preconditioning Jacobian.
   * \param q Q vector of boundary dof.
   * \param qa Qa vector of boundary dof.
   * \param M Preconditioning Jacobian matrix (row order).
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-27
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/StrandSystem/StrandSystem.C lhsPreconJacobian
   */
  virtual void lhsPreconJacobian(const int& npts,
				 const double* q,
				 const double* qa,
				 double* M);

  /**
   * \brief
   * Provides Jacobian of the conserved variables wrt update variables.
   * \param npts Number of points at which to compute Jacobian.
   * \param q Q vector of boundary dof.
   * \param qa Qa vector of boundary dof.
   * \param M Jacobian matrix (row order).
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-11-2
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/StrandSystem/StrandSystem.C lhsTimeJacobian
   */
  virtual void lhsConsVarJacobian(const int& npts,
				  const double* q,
				  const double* qa,
				  double* M);



  /**
   * \brief
   * Provides Jacobian of the primitive variables wrt update variables.
   * \param npts Number of points at which to compute Jacobian.
   * \param q Q vector of boundary dof.
   * \param qa Qa vector of boundary dof.
   * \param M Jacobian matrix (row order).
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-27
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/StrandSystem/StrandSystem.C lhsPrimVarJacobian
   */
  virtual void lhsPrimVarJacobian(const int& npts,
				  const double* q,
				  const double* qa,
				  double* M);

  /**
   * \brief
   * Initializes output at each physical time step.
   * \param step Physical time step.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-25
   * \par Further Documentation:
   */
  virtual void outputInitialize(const int& step)=0;

  /**
   * \brief
   * Initializes output at each physical time step.
   * \param ffile Physical time step.
   * \param npts number of points.
   * \param q Q vector.
   * \param qa Qa vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-25
   * \par Further Documentation:
   */
  virtual void outputSolution(ofstream& ffile,
			      const int& npts,
			      const double* q,
			      const double* qa)=0;

  /**
   * \brief
   * Computes surface forces for surface faces given q on those faces.
   * \param npts Number of points at which to set q.
   * \param tag Tag number of the surface dofs.
   * \param facu Face area components, laid out in (Ax,Ay) pairs.
   * \param q Q vector at the surface faces.
   * \param qa Additional Q vector at the surface faces.
   * \param qx Gradient of Q vector at the surface faces.
   * \param qax Gradient of additional Q vector at the surface faces.
   * \param force Force components at each of the surface faces, laid out in
   * (Fx,Fy) pairs.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-25
   * \par Further Documentation:
   */
  virtual void outputSurfaceForces(const int& npts,
				   const int* tag,
				   const double* A,
				   const double* q,
				   const double* qa,
				   const double* qx,
				   const double* qax,
				   double* force)=0;

  /**
   * \brief
   * Finalizes the StrandSystem class instance.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-04
   * \par Further Documentation:
   */
  virtual void finalize()=0;


 protected:

  int iPrint;
  int iTest;
  int iDebug;
  int nq;
  int nqa;
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
