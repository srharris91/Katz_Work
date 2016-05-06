/**
 * \brief
 * Class SPLamSolution holds the data and specifies the operations to
 * initialize the flow and source terms in various ways.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-25
 */


#ifndef included_SPLamSolution
#define included_SPLamSolution

#include "STRAND_defs.h"
#include "SPLamBc.h"
#include "State.h"
#include "Transport.h"


class SPLamSolution
{

 public:

  /**
   * \brief
   * Constructor for the SPLamSolution class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/SPLamSolution.C SPLamSolution
   */
  SPLamSolution();

  /**
   * \brief
   * Destructor for the SPLamSolution class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/SPLamSolution.C ~SPLamSolution
   */
  ~SPLamSolution();

  /**
   * \brief
   * Initialize the SPLamSolution class and derived classes.
   * \param nq Total number of equations or Q variables in system.
   * \param nqa Total number of additional variables to be stored.
   * \param ndim Number of Spatial Dimensions.
   * \param inviscid Flag to add inviscid terms (0 do not add, 1 add).
   * \param viscous Flag to add viscous terms (0 do not add, 1 add).
   * \param source Flag to add source terms (0 do not add, 1 physical source,
   * -1 MMS source).
   * \param dissipation Flag to add dissipation (0 do not add, 1 add).
   * \param p0 Reference pressure.
   * \param u0 Reference x-velocity.
   * \param v0 Reference y-velocity.
   * \param t0 Reference temperature.
   * \param Re Reynolds number.
   * \param ReRefLength Length used in Reynolds number definition.
   * \param bc pointer to the vector of boundary condition instances.
   * \param state pointer to the state instance.
   * \param transport pointer to the vector of transport instances.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-25
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/SPLamSolution.C initialize
   */
  void initialize(const int& nq,
		  const int& nqa,
		  const int& ndim,
		  const int& inviscid,
		  const int& viscous,
		  const int& source,
		  const int& dissipation,
		  const double& p0,
		  const double& u0,
		  const double& v0,
		  const double& t0,
		  const double& Re,
		  const double& ReRefLength,
		  SPLamBc** bc,
		  State* state,
		  Transport* transport);

  /**
   * \brief
   * Initialize the flow in the derived classes.
   * \param npts Number of points at which to initialize the flow.
   * \param x Dof locations, (x,y) pairs.
   * \param q State variable vector.
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
			double* q) = 0;

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
   * \par Source Code:
   * \snippet src/System/SystemSPLam/SPLamSolution.C initSource
   */
  virtual void initSource(const int& npts,
			  const double* x,
			  double* s);

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
   * \par Source Code:
   * \snippet src/System/SystemSPLam/SPLamSolution.C initSourceBnd
   */
  virtual void initSource(const int& npts,
			  const int* tag,
			  const double* x,
			  double* s);
  
  /**
   * \brief
   * Release memory for the SPLamSolution class and derived classes.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/SPLamSolution.C finalize
   */
  void finalize();


 protected:

  int nq;
  int nqa;
  int ndim;
  int inviscid;
  int viscous;
  int source;
  int dissipation;
  double p0;
  double u0;
  double v0;
  double t0;
  double Re;
  double ReRefLength;
  SPLamBc** bc;
  State* state;
  Transport* transport;


 private:

};
#endif
