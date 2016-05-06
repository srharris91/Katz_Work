/**
 * \brief
 * Class System holds the data and specifies the operations for various
 * systems of equations, such as single phase, turbulent, or multi-species
 * reacting gas.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-11
 */


#ifndef included_System
#define included_System

#include "STRAND_defs.h"
#include "State.h"
#include "Transport.h"
#include "Solution.h"


class System
{
 public:

  /**
   * \brief
   * Constructor for the System class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/System.C System
   */
  System();

  /**
   * \brief
   * Destructor for the System class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/System.C ~System
   */
  ~System();

  /**
   * \brief
   * Reads inputs for the System layer.
   * \param inputFile name of System input file.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/System.C inputRead
   */
  virtual void inputRead(const string& inputFile);

  /**
   * \brief
   * Provides relevant System data to the Numerics routines
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
   * \param nBPatches Number of boundary patches.
   * \param iqgrad Dimension (tmp), flag for whether gradient of Q is
   * required: 0 or 1 for each.
   * \param iqagrad Dimension (tmp), flag for whether gradient of Qa is
   * required: 0 or 1 for each.
   * \param dlim Constant used in limiter computations.
   * \param dqNorm normalization values for Q.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/System.C prepSetup
   */
  virtual void prepSetup(const int& iPrint,
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
   */
  virtual void initSource(const int& npts,
			  const double* x,
			  double* s) = 0;

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
   */
  virtual void initSource(const int& npts,
			  const int* tag,
			  const double* x,
			  double* s) = 0;


 protected:

  State state;
  Transport transport;
  Solution solution;
  int iPrint;
  int iTest;
  int iDebug;
  int nq;
  int nqa;
  int ndim;
  int ncomp;
  int isolution;
  int nBpatches;
  int inviscid;
  int viscous;
  int source;
  int dissipation;


 private:

};
#endif
