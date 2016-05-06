/**
 * \brief
 * Class Strand2dFCSPLamBc holds the data and specifies the operations for
 * various boundary conditions relating to single-phase inviscid or
 * viscous laminar flows.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-12
 */


#ifndef included_Strand2dFCSPLamBc
#define included_Strand2dFCSPLamBc

#include "STRAND2DFC_defs.h"
#include "State.h"
#include "Transport.h"
#include "Solution.h"


class Strand2dFCSPLamBc
{

 public:

  /**
   * \brief
   * Constructor for the Strand2dFCSPLamBc class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/Strand2dFCSPLamBc.C Strand2dFCSPLamBc
   */
  Strand2dFCSPLamBc();

  /**
   * \brief
   * Destructor for the Strand2dFCSPLamBc class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/Strand2dFCSPLamBc.C ~Strand2dFCSPLamBc
   */
  virtual ~Strand2dFCSPLamBc();

  /**
   * \brief
   * Initialize the Strand2dFCSPLamBc class and derived classes.
   * \param nq0 Total number of equations or Q variables in system.
   * \param nqa0 Total number of additional variables to be stored.
   * \param ndim0 Number of spatial dimensions.
   * \param ncomp0 Number of fluid components.
   * \param viscous0 Viscosity flag for viscous flows.
   * \param bValue0 Pointer to array of boundary reference values (p,u,v,T).
   * \param state0 Pointer to the state instance.
   * \param transport0 Pointer to the transport instance.
   * \param solution0 Pointer to the solution instance.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/Strand2dFCSPLamBc.C initialize
   */
  void initialize(const int& nq0,
		  const int& nqa0,
		  const int& ndim0,
		  const int& ncomp0,
		  const int& viscous0,
		  const double& gamma0,
		  const double& rGas0,
		  double* bValue0,
		  State* state0,
		  Transport* transport0,
		  Solution* solution0);

  /**
   * \brief
   * Provides BC vector for a single instance.
   * (Neumann gradient conditions to be added soon)
   * \param nx Outward boundary normal in (nx,ny) pair.
   * \param q Qa vector of boundary dof.
   * \param qa Qa vector of boundary dof.
   * \param rb BC vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-27
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/Strand2dFCSPLamBc.C BCVector
   */
  virtual void BCVector(const double* nx,
			const double* q,
			const double* qa,
			double* rb);

  /**
   * \brief
   * Provides boundary Jacobian contribution for the boundary dof itself.
   * \param nx Outward boundary normal in (nx,ny) pair.
   * \param q Qa vector of boundary dof.
   * \param qa Qa vector of boundary dof.
   * \param M Boundary condition Jacobian matrix (row order) for self.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-30
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/Strand2dFCSPLamBc.C BCVectorSelfJacobian
   */
  virtual void BCVectorSelfJacobian(const double* nx,
				    const double* q,
				    const double* qa,
				    double* M);

  /**
   * \brief
   * Provides boundary Jacobian contribution for the interior dof.
   * \param nx Outward boundary normal in (nx,ny) pair.
   * \param q Qa vector of boundary dof.
   * \param qa Qa vector of boundary dof.
   * \param L Boundary condition Selection matrix.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-30
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/Strand2dFCSPLamBc.C BCSelectionMatrix
   */
  virtual void BCSelectionMatrix(const double* nx,
				 const double* q,
				 const double* qa,
				 double* L);

  /**
   * \brief
   * Provides the BC penalty term.
   * \param inout Flag to determine incoming or outgoing characteristics.
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
  virtual void BCPenalty(const int& inout,
			 const double* A,
			 const double& Pinv0,
			 const double* q,
			 const double* qa,
			 const double* g,
			 const double* uw,
			 double* rb);

  /**
   * \brief
   * Provides the BC penalty term.
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
  virtual void BCPenaltyVis(const double& Pinv0,
			    const double* q,
			    const double* qa,
			    const double* gv,
			    const double* uw,
			    double* rb);

  /**
   * \brief
   * Provides the BC penalty term Jacobian.
   * \param inout Flag to determine incoming or outgoing characteristics.
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
   * \snippet src/System/SystemSPLam/Strand2dFCSPLamBc.C BCPenaltyJacobian
   */
  virtual void BCPenaltyJacobian(const int& inout,
				 const double* A,
				 const double& Pinv0,
				 const double* q,
				 const double* qa,
				 const double* g,
				 const double* uw,
				 double* M);

  /**
   * \brief
   * Initializes penalty data.
   * \param x Coordinates of each dof, laid out in (x,y) pairs.
   * \param f Solution vector at boundary.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-25
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/Strand2dFCSPLamBc.C penaltyData
   */
  virtual void penaltyData(const double* x,
			   double* f);

  /**
   * \brief
   * Returns surface forces.
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
   * 2012-10-10
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/Strand2dFCSPLamBc.C surfaceForces
   */
  virtual void surfaceForces(const double* xs,
			     const double* ys,
			     const double* q,
			     const double* qa,
			     const double* qx,
			     const double* qy,
			     const double* qax,
			     const double* qay,
			     double* force);

  /**
   * \brief
   * Release memory for the Strand2dFCSPLamBc class and derived classes.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/Strand2dFCSPLamBc.C finalize
   */
  void finalize();


 protected:

  int nq;
  int nqa;
  int ncomp;
  int viscous;
  int ndim;
  double gamma;
  double gm1;
  double ggm1;
  double rGas;
  double Prn;
  double* bValue;
  State* state;
  Transport* transport;
  Solution* solution;


 private:

};
#endif
