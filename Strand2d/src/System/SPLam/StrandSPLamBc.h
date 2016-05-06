/**
 * \brief
 * Class StrandSPLamBc holds the data and specifies the operations for
 * various boundary conditions relating to single-phase inviscid or
 * viscous laminar flows.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-12
 */


#ifndef included_StrandSPLamBc
#define included_StrandSPLamBc

#include "STRAND_defs.h"
#include "State.h"
#include "Transport.h"
#include "Solution.h"


class StrandSPLamBc
{

 public:

  /**
   * \brief
   * Constructor for the StrandSPLamBc class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/StrandSPLamBc.C StrandSPLamBc
   */
  StrandSPLamBc();

  /**
   * \brief
   * Destructor for the StrandSPLamBc class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/StrandSPLamBc.C ~StrandSPLamBc
   */
  ~StrandSPLamBc();

  /**
   * \brief
   * Initialize the StrandSPLamBc class and derived classes.
   * \param nq0 Total number of equations or Q variables in system.
   * \param nqa0 Total number of additional variables to be stored.
   * \param ndim0 Number of spatial dimensions.
   * \param ncomp0 Number of fluid components.
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
   * \snippet src/System/SystemSPLam/StrandSPLamBc.C initialize
   */
  void initialize(const int& nq0,
		  const int& nqa0,
		  const int& ndim0,
		  const int& ncomp0,
		  const double& gamma0,
		  const double& rGas0,
		  double* bValue0,
		  Transport* transport0,
		  Solution* solution0);

  /**
   * \brief
   * Provides BC vector for a single instance of  Q and extrapolated values 
   * (Neumann gradient conditions to be added soon)
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
   * \snippet src/System/SystemSPLam/StrandSPLamBc.C BCVector
   */
  virtual void BCVector(const double* nx,
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
   * \snippet src/System/SystemSPLam/StrandSPLamBc.C BCVectorSelfJacobian
   */
  virtual void BCVectorSelfJacobian(const double* nx,
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
   * \snippet src/System/SystemSPLam/StrandSPLamBc.C BCVectorInteriorJacobian
   */
  virtual void BCVectorInteriorJacobian(const double* nx,
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
   * \snippet src/System/SystemSPLam/StrandSPLamBc.C surfaceForces
   */
  virtual void surfaceForces(const double* A,
			     const double* q,
			     const double* qa,
			     const double* qx,
			     const double* qax,
			     double* force);

  /**
   * \brief
   * Release memory for the StrandSPLamBc class and derived classes.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/StrandSPLamBc.C finalize
   */
  void finalize();


 protected:

  int nq;
  int nqa;
  int ncomp;
  int ndim;
  double gamma;
  double gm1;
  double ggm1;
  double rGas;
  double* bValue;
  Transport* transport;
  Solution* solution;


 private:

};
#endif
