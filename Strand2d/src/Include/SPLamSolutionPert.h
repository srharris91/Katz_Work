/**
 * \brief
 * Class SPLamSolutionPert holds the data and specifies the operations for
 * free stream solution initialization.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-25
 */


#ifndef included_SPLamSolutionPert
#define included_SPLamSolutionPert

#include "STRAND_defs.h"
#include "SPLamSolution.h"


class SPLamSolutionPert: public SPLamSolution
{

 public:

  /**
   * \brief
   * Constructor for the SPLamSolutionPert class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/SPLamSolutionPert.C SPLamSolutionPert
   */
  SPLamSolutionPert();

  /**
   * \brief
   * Destructor for the SPLamSolutionPert class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/SPLamSolutionPert.C ~SPLamSolutionPert
   */
  ~SPLamSolutionPert();

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
   * \snippet src/System/SystemSPLam/SPLamSolutionPert.C initFlow
   */
  void initFlow(const int& npts,
		const double* x,
		double* q);


 private:

};
#endif
