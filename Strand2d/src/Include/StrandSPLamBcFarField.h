/**
 * \brief
 * Class StrandSPLamBcFarField holds the data and specifies the operations for
 * a far field boundary condition for single-phase inviscid or
 * viscous laminar flows.
 *
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-12
 */


#ifndef included_StrandSPLamBcFarField
#define included_StrandSPLamBcFarField

#include "STRAND_defs.h"
#include "StrandSPLamBc.h"


class StrandSPLamBcFarField: public StrandSPLamBc
{

 public:

  /**
   * \brief
   * Constructor for the StrandSPLamBcFarField class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/StrandSPLamBcFarField.C StrandSPLamBcFarField
   */
  StrandSPLamBcFarField();

  /**
   * \brief
   * Destructor for the StrandSPLamBcFarField class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/StrandSPLamBcFarField.C ~StrandSPLamBcFarField
   */
  ~StrandSPLamBcFarField();

 private:

};
#endif
