/**
 * \brief
 * Class Perturb holds the data and specifies the operations for
 * solution initialization using perturbed free stream conditions.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-13
 */


#ifndef included_Perturb
#define included_Perturb

#include "PHYSICS_defs.h"
#include "SolutionType.h"


class Perturb: public SolutionType
{

 public:

  /**
   * \brief
   * Constructor for the Perturb class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Perturb.C Perturb
   */
  Perturb();

  /**
   * \brief
   * Destructor for the Perturb class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-13
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Perturb.C ~Perturb
   */
  ~Perturb();

  /**
   * \brief
   * Reads inputs and initializes the Perturb class objects.
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
   * \snippet src/Solution/Perturb.C initialize
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
   * \snippet src/Solution/Perturb.C getQ
   */
  void getQ(const int& npts,
	    const double* x,
	    double* q);

  /**
   * \brief
   * Deallocates any memory for Perturb class objects.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-27
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/Solution/Perturb.C finalize
   */
  void finalize();


 private:

  double pert;
};
#endif
