/**
 * \brief
 * Class SPLamSolutionMMS holds the data and specifies the operations for
 * manufactured solution initialization.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-25
 */


#ifndef included_SPLamSolutionMMS
#define included_SPLamSolutionMMS

#include "STRAND_defs.h"
#include "SPLamSolution.h"


class SPLamSolutionMMS: public SPLamSolution
{

 public:

  /**
   * \brief
   * Constructor for the SPLamSolutionMMS class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/SPLamSolutionMMS.C SPLamSolutionMMS
   */
  SPLamSolutionMMS();

  /**
   * \brief
   * Destructor for the SPLamSolutionMMS class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/SPLamSolutionMMS.C ~SPLamSolutionMMS
   */
  ~SPLamSolutionMMS();

  /**
   * \brief
   * Initializes free stream flow.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-25
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/SPLamSolutionMMS.C initFlow
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
   * \snippet src/System/SystemSPLam/SPLamSolutionMMS.C initSource
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
   * \snippet src/System/SystemSPLam/SPLamSolutionMMS.C initSourceBnd
   */
  void initSource(const int& npts,
		  const int* tag,
		  const double* x,
		  double* s);


 private:
  /*
  const double c1    = 100000.;
  const double cx1   = 20000.;
  const double cy1   = 17500.;
  const double cxy1  =-25000.;
  const double ax1   = .33;
  const double ay1   = .42;
  const double axy1  = .0363;
  const double c2    = 70.;
  const double cx2   = 7.;
  const double cy2   =-8.;
  const double cxy2  = 5.5;
  const double ax2   = .45;
  const double ay2   = .48;
  const double axy2  = .0069;
  const double c3    = 90.;
  const double cx3   =-5.;
  const double cy3   = 10.;
  const double cxy3  =-11.;
  const double ax3   = .45;
  const double ay3   = .46;
  const double axy3  = .0294;
  const double c4    = 300.;
  const double cx4   = 20.;
  const double cy4   = 17.;
  const double cxy4  =-25.;
  const double ax4   = .33;
  const double ay4   = .42;
  const double axy4  = .0363;
  */
};
#endif
