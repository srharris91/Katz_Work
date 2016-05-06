/**
 * \brief
 * Class Solution is responsible for initializing the solution based on
 * various options, inclucing manufactured solutions, exact solutions,
 * or file inputs.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-27
 */


#ifndef included_Solution
#define included_Solution

#include "PHYSICS_defs.h"
#include "State.h"
#include "Transport.h"
#include "SolutionType.h"
#include "FreeStream.h"
#include "Perturb.h"
#include "Mms2d.h"
#include "Mms3d.h"


class Solution
{
 public:

  /**
   * \brief
   * Constructor for the Solution class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Solution.C Solution
   */
  Solution();

  /**
   * \brief
   * Destructor for the Solution class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-27
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Solution.C ~Solution
   */
  ~Solution();

  /**
   * \brief
   * Reads inputs and initializes the Solution class objects.
   * \param isolution Flag to determine solution type.
   * \param nq Dimension of q (number of equations).
   * \param ndim Number of spatial dimensions.
   * \param inputFile Name of Solution input file.
   * \param state Pointer to the current State instance.
   * \param transport Pointer to the current Transport instance.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-28
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Solution.C initialize
   */
  void initialize(const int& isolution,
		  const int& nq,
		  const int& ndim,
		  const string& inputFile,
		  State* state,
		  Transport* transport);
  void initialize(const int& isolution,
		  const int& nq,
		  const int& ndim,
		  const string& inputFile);

  /**
   * \brief
   * Returns a pointer to the array of solution reference values.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-02
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Solution.C initialize
   */
  double* getRefValues();

  /**
   * \brief
   * Returns the q vector given the location in space.
   * \param npts Number of points at which to return q.
   * \param x Coordinates of q.
   * \param q q-vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-02
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Solution.C getQ
   */
  void getQ(const int& npts,
	    const double* x,
	    double* q);

  /**
   * \brief
   * Returns the x-derivative of the q vector given the location in space.
   * \param npts Number of points at which to return qx.
   * \param x Coordinates of qx.
   * \param qx qx-vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-02
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Solution.C getQx
   */
  void getQx(const int& npts,
	     const double* x,
	     double* qx);

  /**
   * \brief
   * Returns the y-derivative of the q vector given the location in space.
   * \param npts Number of points at which to return qy.
   * \param x Coordinates of qy.
   * \param qy qy-vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-02
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Solution.C getQy
   */
  void getQy(const int& npts,
	     const double* x,
	     double* qy);

  /**
   * \brief
   * Returns the z-derivative of the q vector given the location in space.
   * \param npts Number of points at which to return qz.
   * \param x Coordinates of qz.
   * \param qz qz-vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-02
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Solution.C getQz
   */
  void getQz(const int& npts,
	     const double* x,
	     double* qz);

  /**
   * \brief
   * Returns the second x-derivative of the q vector given the location
   * in space.
   * \param npts Number of points at which to return qxx.
   * \param x Coordinates of qxx.
   * \param qxx qxx-vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-02
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Solution.C getQxx
   */
  void getQxx(const int& npts,
	      const double* x,
	      double* qxx);

  /**
   * \brief
   * Returns the second y-derivative of the q vector given the location
   * in space.
   * \param npts Number of points at which to return qyy.
   * \param x Coordinates of qyy.
   * \param qyy qyy-vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-02
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Solution.C getQyy
   */
  void getQyy(const int& npts,
	      const double* x,
	      double* qyy);

  /**
   * \brief
   * Returns the second z-derivative of the q vector given the location
   * in space.
   * \param npts Number of points at which to return qzz.
   * \param x Coordinates of qzz.
   * \param qzz qzz-vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-02
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Solution.C getQzz
   */
  void getQzz(const int& npts,
	      const double* x,
	      double* qzz);

  /**
   * \brief
   * Returns the mixed xy-derivative of the q vector given the location
   * in space.
   * \param npts Number of points at which to return qxy.
   * \param x Coordinates of qxy.
   * \param qxy qxy-vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-02
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Solution.C getQxy
   */
  void getQxy(const int& npts,
	      const double* x,
	      double* qxy);

  /**
   * \brief
   * Returns the mixed xz-derivative of the q vector given the location
   * in space.
   * \param npts Number of points at which to return qxz.
   * \param x Coordinates of qxz.
   * \param qxz qxz-vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-02
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Solution.C getQxz
   */
  void getQxz(const int& npts,
	      const double* x,
	      double* qxz);

  /**
   * \brief
   * Returns the mixed yz-derivative of the q vector given the location
   * in space.
   * \param npts Number of points at which to return qyz.
   * \param x Coordinates of qyz.
   * \param qyz qyz-vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-02
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Solution.C getQyz
   */
  void getQyz(const int& npts,
	      const double* x,
	      double* qyz);

  /**
   * \brief
   * Deallocates any memory used.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-27
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Solution.C finalize
   */
  void finalize();


 protected:

  SolutionType* solutionType;


 private:

};
#endif
