/**
 * \brief
 * Class SPLamBc holds the data and specifies the operations for
 * various boundary conditions relating to single-phase inviscid or
 * viscous laminar flows.
 * \author
 * Aaron Katz
 * \version
 * 1.0
 * \date
 * 2012-09-12
 */


#ifndef included_SPLamBc
#define included_SPLamBc

#include "STRAND_defs.h"
#include "State.h"
#include "Transport.h"
#include "Solution.h"


class SPLamBc
{

 public:

  /**
   * \brief
   * Constructor for the SPLamBc class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/SPLamBc.C SPLamBc
   */
  SPLamBc();

  /**
   * \brief
   * Destructor for the SPLamBc class.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-12
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/SPLamBc.C ~SPLamBc
   */
  ~SPLamBc();

  /**
   * \brief
   * Initialize the SPLamBc class and derived classes.
   * \param nq0 Total number of equations or Q variables in system.
   * \param nqa0 Total number of additional variables to be stored.
   * \param ncomp0 Number of fluid components.
   * \param bValue0 Pointer to array of boundary reference values (p,u,v,T).
   * \param state0 Pointer to the state instance.
   * \param transport0 Pointer to the transport instance.
   * \param solution0 Pointer to the solution instance.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   */
  void initialize(const int& nq0,
		  const int& nqa0,
		  const int& ncomp0,
		  double* bValue0,
		  State* state0,
		  Transport* transport0,
		  Solution* solution0);

  /**
   * \brief
   * Release memory for the SPLamBc class and derived classes.
   * \author
   * Aaron Katz
   * \version
   * 1.0
   * \date
   * 2012-09-11
   * \par Further Documentation:
   * \par Source Code:
   * \snippet src/System/SystemSPLam/SPLamBc.C finalize
   */
  void finalize();


 protected:

  int nq;
  int nqa;
  int ncomp;
  double* bValue;
  State* state;
  Transport* transport;
  Solution* solution;


 private:

};
#endif
