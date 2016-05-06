/**
 * \brief
 * Class SPLamBcInflow holds the data and specifies the operations for
 * an inflow boundary condition for single-phase inviscid or
 * viscous laminar flows.
 *
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-12
 */


#ifndef included_SPLamBcInflow
#define included_SPLamBcInflow

#include "STRAND_defs.h"
#include "SPLamBc.h"


class SPLamBcInflow: public SPLamBc
{

 public:

  /**
   * \brief
   * Constructor for the SPLamBcInflow class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/SPLamBcInflow.C SPLamBcInflow
   */
  SPLamBcInflow();

  /**
   * \brief
   * Destructor for the SPLamBcInflow class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/SPLamBcInflow.C ~SPLamBcInflow
   */
  ~SPLamBcInflow();

 private:

};
#endif
