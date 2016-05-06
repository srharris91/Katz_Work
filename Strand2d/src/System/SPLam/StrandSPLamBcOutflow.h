/**
 * \brief
 * Class StrandSPLamBcOutflow holds the data and specifies the operations for
 * an outflow boundary condition for single-phase inviscid or
 * viscous laminar flows.
 *
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-12
 */


#ifndef included_StrandSPLamBcOutflow
#define included_StrandSPLamBcOutflow

#include "STRAND_defs.h"
#include "StrandSPLamBc.h"


class StrandSPLamBcOutflow: public StrandSPLamBc
{

 public:

  /**
   * \brief
   * Constructor for the StrandSPLamBcOutflow class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/StrandSPLamBcOutflow.C StrandSPLamBcOutflow
   */
  StrandSPLamBcOutflow();

  /**
   * \brief
   * Destructor for the StrandSPLamBcOutflow class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/StrandSPLamBcOutflow.C ~StrandSPLamBcOutflow
   */
  ~StrandSPLamBcOutflow();

  /**
   * \brief
   * Provides outflow BC vector for a single instance of Q and
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
   * src/System/SystemSPLam/StrandSPLamBcInflow.C BCVector
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
   * \snippet src/System/SystemSPLam/StrandSPLamBcOutflow.C BCVectorSelfJacobian
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
   * \snippet src/System/SystemSPLam/StrandSPLamBcOutflow.C BCVectorInteriorJacobian
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
