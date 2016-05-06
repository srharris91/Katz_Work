/**
 * \brief
 * Implementation file for index functions for multidimensional arrays used
 * StrandBlockSolver.
 *
 * \author
 * Aaron Katz
 *
 * \version
 * 1.0
 *
 * \date
 * 2012-09-20
 */


#include "StrandBlockSolver.h"

void StrandBlockSolver::indlsp(const int& k, const int& j, const int& n,
			       int& i){
  i = 2*((nPstr+1)*n+j)+k;}
