/**
 * \brief
 * Class SystemSPLam
 * holds the data and specifies the operations for the single phase inviscid
 * or laminar viscous system of equations.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-11
 */

#ifndef included_SystemSPLam
#define included_SystemSPLam

#include "System.h"
#include "SPLamBc.h"
#include "SPLamBcInviscidWall.h"
#include "SPLamBcViscousWall.h"
#include "SPLamBcInflow.h"
#include "SPLamBcOutflow.h"
#include "SPLamBcFarField.h"
#include "SPLamBcDirichlet.h"
#include "SPLamBcFrozen.h"


class SystemSPLam: public System
{
 public:

  /**
   * \brief
   * Constructor for the SystemSPLam class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/SystemSPLam.C SystemSPLam
   */
  SystemSPLam();

  /**
   * \brief
   * Destructor for the SystemSPLam class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/SystemSPLam.C ~SystemSPLam
   */
  ~SystemSPLam();

  /**
   * \brief
   * Reads inputs for the System layer, and allocates instances of
   * boundary conditions and Physics classes.
   * \param inputFile name of System input file.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \include src/System/SystemSPLam/inputRead.C
   */
  void inputRead(const string& inputFile);

  /**
   * \brief
   * Provides relevant Aystem data to the Numerics routines
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
   * \param source Flag to add source terms (0 do not add, 1 physical source,
   * -1 MMS source).
   * \param dissipation Flag to add dissipation (0 do not add, 1 add).
   * \param nBpatches Number of boundary patches.
   * \param iqgrad Dimension (tmp), flag for whether gradient of Q is
   * required: 0 or 1 for each.
   * \param iqagrad Dimension (tmp), flag for whether gradient of Qa is
   * required: 0 or 1 for each.
   * \param dlim Constant used in limiter computations.
   * \param rmsNorm normalization values for Q.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \include src/System/SystemSPLam/prepSetup.C
   */
  void prepSetup(const int& iPrint,
		 const int& iTest,
		 const int& iDebug,
		 const int& tmp,
		 int& nq,
		 int& nqa,
		 int& ndim,
		 int& inviscid,
		 int& viscous,
		 int& source,
		 int& dissipation,
		 int& nBpatches,
		 int* iqgradT,
		 int* iqagradT,
		 double& dlim,
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
   * \par Source Code:
   * \include src/System/SystemSPLam/initFlow.C
   */
  void initFlow(const int& npts,
		const double* x,
		double* q);

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
   * \include src/System/SystemSPLam/initSource.C
   */
  void initSource(const int& npts,
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
   * \include src/System/SystemSPLam/initSource.C
   */
  void initSource(const int& npts,
		  const int* tag,
		  const double* x,
		  double* s);


 private:

  SPLamBc** bc;
};
#endif
