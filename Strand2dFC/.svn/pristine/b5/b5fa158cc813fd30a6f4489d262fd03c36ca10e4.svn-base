/**
 * \brief
 * Class Strand2dFCSPTurbSABcViscousWall holds the data and specifies the
 * operations for
 * a dirichlet boundary condition for single-phase inviscid or
 * viscous laminar flows.
 *
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2013-01-02
 */


#ifndef included_Strand2dFCSPTurbSABcViscousWall
#define included_Strand2dFCSPTurbSABcViscousWall

#include "STRAND2DFC_defs.h"
#include "Strand2dFCSPTurbSABc.h"


class Strand2dFCSPTurbSABcViscousWall: public Strand2dFCSPTurbSABc
{

 public:

  /**
   * \brief
   * Constructor for the Strand2dFCSPTurbSABcViscousWall class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-01-02
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPTurbSA/Strand2dFCSPTurbSABcViscousWall.C Strand2dFCSPTurbSABcViscousWall
   */
  Strand2dFCSPTurbSABcViscousWall();

  /**
   * \brief
   * Destructor for the Strand2dFCSPTurbSABcViscousWall class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-01-02
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPTurbSA/Strand2dFCSPTurbSABcViscousWall.C ~Strand2dFCSPTurbSABcViscousWall
   */
  ~Strand2dFCSPTurbSABcViscousWall();

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
   * \snippet src/System/SystemSPTurbSA/Strand2dFCSPTurbSABcViscousWall.C
   */
  void BCVector(const double* nx,
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
   * \snippet src/System/SystemSPTurbSA/Strand2dFCSPTurbSABcViscousWall.C
   */
  void BCVectorSelfJacobian(const double* nx,
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
   * \snippet src/System/SystemSPTurbSA/Strand2dFCSPTurbSABcViscousWall.C
   */
  void BCSelectionMatrix(const double* nx,
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
   * \snippet src/System/SystemSPTurbSA/Strand2dFCSPTurbSABcViscousWall.C
   */
  void BCPenalty(const int& inout,
		 const double* A,
		 const double& Pinv0,
		 const double* q,
		 const double* qa,
		 const double* g,
		 const double* gv,
		 const double* uw,
		 double* rb);


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
   * \snippet src/System/SystemSPTurbSA/Strand2dFCSPTurbSABcViscousWall.C surfaceForces
   */
  void surfaceForces(const double* xs,
		     const double* ys,
		     const double* q,
		     const double* qa,
		     const double* qx,
		     const double* qy,
		     const double* qax,
		     const double* qay,
		     double* force);


 private:

};
#endif
