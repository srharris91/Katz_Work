/**
 * \brief
 * Class SPLamSolutionFS holds the data and specifies the operations for
 * free stream solution initialization.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-25
 */


#ifndef included_SPLamSolutionFS
#define included_SPLamSolutionFS

#include "STRAND_defs.h"
#include "SPLamSolution.h"


class SPLamSolutionFS: public SPLamSolution
{

 public:

  /**
   * \brief
   * Constructor for the SPLamSolutionFS class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/SPLamSolutionFS.C SPLamSolutionFS
   */
  SPLamSolutionFS();

  /**
   * \brief
   * Destructor for the SPLamSolutionFS class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/SPLamSolutionFS.C ~SPLamSolutionFS
   */
  ~SPLamSolutionFS();

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
   * \snippet src/System/SystemSPLam/SPLamSolutionFS.C initFlow
   */
  void initFlow(const int& npts,
		const double* x,
		double* q);


 private:

};
#endif
