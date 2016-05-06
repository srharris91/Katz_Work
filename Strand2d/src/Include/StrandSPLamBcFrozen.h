/**
 * \brief
 * Class StrandSPLamBcFrozen holds the data and specifies the operations for
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


#ifndef included_StrandSPLamBcFrozen
#define included_StrandSPLamBcFrozen

#include "STRAND_defs.h"
#include "StrandSPLamBc.h"


class StrandSPLamBcFrozen: public StrandSPLamBc
{

 public:

  /**
   * \brief
   * Constructor for the StrandSPLamBcFrozen class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/StrandSPLamBcFrozen.C StrandSPLamBcFrozen
   */
  StrandSPLamBcFrozen();

  /**
   * \brief
   * Destructor for the StrandSPLamBcFrozen class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/StrandSPLamBcFrozen.C ~StrandSPLamBcFrozen
   */
  ~StrandSPLamBcFrozen();

 private:

};
#endif
