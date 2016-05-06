/**
 * \brief
 * Class StrandSPTurbSABcDirichlet holds the data and specifies the operations for
 * a dirichlet boundary condition for single-phase inviscid or
 * viscous Turbulent flows.
 *
 * \author
 * Aaron Katz
 * \version
 * 1.0 T
 * \date
 * 2012-11-12
 */


#ifndef included_StrandSPTurbSABcDirichlet
#define included_StrandSPTurbSABcDirichlet

#include "STRAND_defs.h"
#include "StrandSPTurbSABc.h"


class StrandSPTurbSABcDirichlet: public StrandSPTurbSABc
{

 public:

  /**
   * \brief
   * Constructor for the StrandSPTurbSABcDirichlet class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPTurbSA/StrandSPTurbSABcDirichlet.C StrandSPTurbSABcDirichlet
   */
  StrandSPTurbSABcDirichlet();

  /**
   * \brief
   * Destructor for the StrandSPTurbSABcDirichlet class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPTurbSA/StrandSPTurbSABcDirichlet.C ~StrandSPTurbSABcDirichlet
   */
  ~StrandSPTurbSABcDirichlet();

  /**
   * \brief
   * Provides Dirichlet BC vector for a single instance of Q and
   * extrapolated values.
   * \param nx Outward boundary normal in (nx,ny) pair.
   * \param wx Wall velocity in (wx,wy) pair.
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
   * \par Source Code:
   * \snippet
   * src/System/SystemSPTurbSA/StrandSPTurbSABcDirichlet.C BCVector
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
   * \snippet src/System/SystemSPTurbSA/StrandSPTurbSABcDirichlet.C BCVectorSelfJacobian
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
   * \snippet src/System/SystemSPTurbSA/StrandSPTurbSABcDirichlet.C BCVectorInteriorJacobian
   */
  void BCVectorInteriorJacobian(const double* nx,
				const double* qe,
				const double* qae,
				const double* q,
				const double* qa,
				double* M);


 private:

};
#endif
