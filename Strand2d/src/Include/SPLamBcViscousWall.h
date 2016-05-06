/**
 * \brief
 * Class SPLamBcViscousWall holds the data and specifies the operations for
 * a viscous wall boundary condition for single-phase inviscid or
 * viscous laminar flows.
 *
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-12
 */


#ifndef included_SPLamBcViscousWall
#define included_SPLamBcViscousWall

#include "STRAND_defs.h"
#include "SPLamBc.h"


class SPLamBcViscousWall: public SPLamBc
{

 public:

  /**
   * \brief
   * Constructor for the SPLamBcViscousWall class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/SPLamBcViscousWall.C SPLamBcViscousWall
   */
  SPLamBcViscousWall();

  /**
   * \brief
   * Destructor for the SPLamBcViscousWall class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/SPLamBcViscousWall.C ~SPLamBcViscousWall
   */
  ~SPLamBcViscousWall();

 private:

};
#endif
