/**
 * \brief
 * Class Mms2d holds the data and specifies the operations for
 * solution initialization using a 2d trigonometric manufactured solution.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-13
 */


#ifndef included_Mms2d
#define included_Mms2d

#include "PHYSICS_defs.h"
#include "SolutionType.h"


class Mms2d: public SolutionType
{

 public:

  /**
   * \brief
   * Constructor for the Mms2d class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Mms2d.C Mms2d
   */
  Mms2d();

  /**
   * \brief
   * Destructor for the Mms2d class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Mms2d.C ~Mms2d
   */
  ~Mms2d();

  /**
   * \brief
   * Reads inputs and initializes the Mms2d class objects.
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
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Mms2d.C initialize
   */
  void initialize(const int& nq,
		  const int& ndim,
		  const string& inputFile,
		  State* state,
		  Transport* transport);
  void initialize(const int& nq,
		  const int& ndim,
		  const string& inputFile);

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
   * \snippet src/Solution/Mms2d.C getQ
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
   * \snippet src/Solution/Mms2d.C getQx
   */
  void getQx(const int& npts,
	     const double* x,
	     double* qx);

  /**
   * \brief
   * Returns the y-derivative of the q vector given the location in space.
   * \param npts Number of points at which to return qy.
   * \param y Coordinates of qy.
   * \param qy qy-vector.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-10-02
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Mms2d.C getQy
   */
  void getQy(const int& npts,
	     const double* x,
	     double* qy);

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
   * \snippet src/Solution/Mms2d.C getQxx
   */
  void getQxx(const int& npts,
	      const double* x,
	      double* qxx);

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
   * \snippet src/Solution/Mms2d.C getQxy
   */
  void getQxy(const int& npts,
	      const double* x,
	      double* qxy);

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
   * \snippet src/Solution/Mms2d.C getQyy
   */
  void getQyy(const int& npts,
	      const double* x,
	      double* qyy);

  /**
   * \brief
   * Deallocates any memory for Mms2d class objects.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-27
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Mms2d.C finalize
   */
  void finalize();


 private:

  double period;
  double amplitude;
  double *ax,*ay,*axy,*bx,*by,*bxy,*cx,*cy,*cxy;
};
#endif
