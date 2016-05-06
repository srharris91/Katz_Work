#include "Strand2dFCSPTurbSA.h"
  
  
void Strand2dFCSPTurbSA::initWallDist(const int& npts,
                                      const double* dw,
                                      double* qa)
 {
   int iqa;
   for (int n=0; n<npts; n++){
     iqa       = nqa*n;
     qa[iqa+7] = dw[n];
   }
 }
