/**
 * \brief
 * Class StrandSPTurbSABcInviscidWall holds the data and specifies the
 * operations for
 * an inviscid wall boundary condition for single-phase inviscid or
 * viscous turbulent flows.
 *
 * \author
 * Aaron Katz
 * \version
 * 1.0 T
 * \date
 * 2012-11-12
 */


#ifndef included_StrandSPTurbSABcInviscidWall
#define included_StrandSPTurbSABcInviscidWall

#include "STRAND_defs.h"
#include "StrandSPTurbSABc.h"


class StrandSPTurbSABcInviscidWall: public StrandSPTurbSABc
{

 public:

  /**
   * \brief
   * Constructor for the StrandSPTurbSABcInviscidWall class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPTurbSA/StrandSPTurbSABcInviscidWall.C StrandSPTurbSABcInviscidWall
   */
  StrandSPTurbSABcInviscidWall();

  /**
   * \brief
   * Destructor for the StrandSPTurbSABcInviscidWall class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPTurbSA/StrandSPTurbSABcInviscidWall.C ~StrandSPTurbSABcInviscidWall
   */
  ~StrandSPTurbSABcInviscidWall();

  /**
   * \brief
   * Provides inviscid wall BC vector for a single instance of Q and
   * extrapolated values.
   * \param nx Outward boundary normal in (nx,ny) pair.
   * \param wx Wall velocity in (wx,wy) pair.
   * \param qe Extrapolated Q vector at boundary location.
   * \param qae Extrapolated Qa vector at boundary location.
   * \param q Qa vector of boundary dof.
   * \param qa Qa vector of boundary dof.
   * \param r BC vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-27
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPTurbSA/StrandSPTurbSABcInviscidWall.C BCVector
   */
  void BCVector(const double* nx,
		const double* wx,
		const double* qe,
		const double* qae,
		const double* q,
		const double* qa,
		double* r);

  /**
   * \brief
   * Provides boundary Jacobian contribution for the boundary dof itself.
   * \param nx Outward boundary normal in (nx,ny) pair.
   * \param qe Extrapolated Q vector at boundary location.
   * \param qae Extrapolated Qa vector at boundary location.
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
   * \snippet src/System/SystemSPTurbSA/StrandSPTurbSABcInviscidWall.C BCVectorSelfJacobian
   */
  void BCVectorSelfJacobian(const double* nx,
			    const double* qe,
			    const double* qae,
			    const double* q,
			    const double* qa,
			    double* M);

  /**
   * \brief
   * Provides boundary Jacobian contribution for the interior dof.
   * \param nx Outward boundary normal in (nx,ny) pair.
   * \param qe Extrapolated Q vector at boundary location.
   * \param qae Extrapolated Qa vector at boundary location.
   * \param q Qa vector of boundary dof.
   * \param qa Qa vector of boundary dof.
   * \param M Boundary condition Jacobian matrix (row order) for interior.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-30
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPTurbSA/StrandSPTurbSABcInviscidWall.C BCVectorInteriorJacobian
   */
  void BCVectorInteriorJacobian(const double* nx,
				const double* qe,
				const double* qae,
				const double* q,
				const double* qa,
				double* M);

  /**
   * \brief
   * Computes surface forces for a surface face given q on this face.
   * \param facu Face area components, (Ax,Ay) pair.
   * \param q Q vector at the surface face.
   * \param qa Additional Q vector at the surface face.
   * \param qx Gradient of Q vector at the surface face.
   * \param qax Gradient of additional Q vector at the surface face.
   * \param force Force components at the surface face, laid out as
   * (Fx,Fy) pair.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-10
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPTurbSA/StrandSPTurbSABcInviscidWall.C surfaceForces
   */
  void surfaceForces(const double* A,
		     const double* q,
		     const double* qa,
		     const double* qx,
		     const double* qax,
		     double* force);


 private:

};
#endif
