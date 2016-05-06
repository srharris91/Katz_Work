/**
 * \brief
 * Class Strand2dFCSPLamBcOutflow holds the data and specifies the operations for
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


#ifndef included_Strand2dFCSPLamBcOutflow
#define included_Strand2dFCSPLamBcOutflow

#include "STRAND2DFC_defs.h"
#include "Strand2dFCSPLamBc.h"


class Strand2dFCSPLamBcOutflow: public Strand2dFCSPLamBc
{

 public:

  /**
   * \brief
   * Constructor for the Strand2dFCSPLamBcOutflow class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-01-02
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/Strand2dFCSPLamBcOutflow.C Strand2dFCSPLamBcOutflow
   */
  Strand2dFCSPLamBcOutflow();

  /**
   * \brief
   * Destructor for the Strand2dFCSPLamBcOutflow class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-01-02
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/Strand2dFCSPLamBcOutflow.C ~Strand2dFCSPLamBcOutflow
   */
  ~Strand2dFCSPLamBcOutflow();


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
   * \snippet src/System/SystemSPLam/Strand2dFCSPLamBcOutflow.C
   */
  void BCPenalty(const int& inout,
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
   * \snippet src/System/SystemSPLam/Strand2dFCSPLamBcOutflow.C
   */
  void BCPenaltyVis(const double& Pinv0,
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
   * \snippet src/System/SystemSPLam/Strand2dFCSPLamBcOutflow.C
   */
  void BCPenaltyJacobian(const int& inout,
			 const double* A,
			 const double& Pinv0,
			 const double* q,
			 const double* qa,
			 const double* g,
			 const double* uw,
			 double* M);


 private:

};
#endif
