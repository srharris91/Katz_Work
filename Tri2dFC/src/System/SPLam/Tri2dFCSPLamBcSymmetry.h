/**
 * \brief
 * Class Tri2dFCSPLamBcSymmetry holds the data and specifies the operations for
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


#ifndef included_Tri2dFCSPLamBcSymmetry
#define included_Tri2dFCSPLamBcSymmetry

#include "TRI2DFC_defs.h"
#include "Tri2dFCSPLamBc.h"


class Tri2dFCSPLamBcSymmetry: public Tri2dFCSPLamBc
{

 public:

  /**
   * \brief
   * Constructor for the Tri2dFCSPLamBcSymmetry class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-01-02
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/Tri2dFCSPLamBcSymmetry.C Tri2dFCSPLamBcSymmetry
   */
  Tri2dFCSPLamBcSymmetry();

  /**
   * \brief
   * Destructor for the Tri2dFCSPLamBcSymmetry class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2013-01-02
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/Tri2dFCSPLamBcSymmetry.C ~Tri2dFCSPLamBcSymmetry
   */
  ~Tri2dFCSPLamBcSymmetry();

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
   * \snippet src/System/SystemSPLam/Tri2dFCSPLamBcSymmetry.C
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
   * \snippet src/System/SystemSPLam/Tri2dFCSPLamBcSymmetry.C
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
   * \snippet src/System/SystemSPLam/Tri2dFCSPLamBcSymmetry.C
   */
  void BCSelectionMatrix(const double* nx,
			 const double* q,
			 const double* qa,
			 double* L);


 private:

};
#endif
