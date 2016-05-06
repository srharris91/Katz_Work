/**
 * \brief
 * Class SPLamBcFarField holds the data and specifies the operations for
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


#ifndef included_SPLamBcFarField
#define included_SPLamBcFarField

#include "STRAND_defs.h"
#include "SPLamBc.h"


class SPLamBcFarField: public SPLamBc
{

 public:

  /**
   * \brief
   * Constructor for the SPLamBcFarField class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/SPLamBcFarField.C SPLamBcFarField
   */
  SPLamBcFarField();

  /**
   * \brief
   * Destructor for the SPLamBcFarField class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/SPLamBcFarField.C ~SPLamBcFarField
   */
  ~SPLamBcFarField();

 private:

};
#endif
