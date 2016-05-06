/**
 * \brief
 * Class SPLamBcFrozen holds the data and specifies the operations for
 * a frozen boundary condition for single-phase inviscid or
 * viscous laminar flows.
 *
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-12
 */


#ifndef included_SPLamBcFrozen
#define included_SPLamBcFrozen

#include "STRAND_defs.h"
#include "SPLamBc.h"


class SPLamBcFrozen: public SPLamBc
{

 public:

  /**
   * \brief
   * Constructor for the SPLamBcFrozen class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/SPLamBcFrozen.C SPLamBcFrozen
   */
  SPLamBcFrozen();

  /**
   * \brief
   * Destructor for the SPLamBcFrozen class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/SPLamBcFrozen.C ~SPLamBcFrozen
   */
  ~SPLamBcFrozen();

 private:

};
#endif
