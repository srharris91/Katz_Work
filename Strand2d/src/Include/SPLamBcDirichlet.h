/**
 * \brief
 * Class SPLamBcDirichlet holds the data and specifies the operations for
 * a dirichlet boundary condition for single-phase inviscid or
 * viscous laminar flows.
 *
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-12
 */


#ifndef included_SPLamBcDirichlet
#define included_SPLamBcDirichlet

#include "STRAND_defs.h"
#include "SPLamBc.h"


class SPLamBcDirichlet: public SPLamBc
{

 public:

  /**
   * \brief
   * Constructor for the SPLamBcDirichlet class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/SPLamBcDirichlet.C SPLamBcDirichlet
   */
  SPLamBcDirichlet();

  /**
   * \brief
   * Destructor for the SPLamBcDirichlet class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet
   * src/System/SystemSPLam/SPLamBcDirichlet.C ~SPLamBcDirichlet
   */
  ~SPLamBcDirichlet();

 private:

};
#endif
