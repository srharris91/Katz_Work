#include "Strand2dFCSPLam.h"


void Strand2dFCSPLam::bcConflict(const int& npts,
			      int* tag,
			      const int* tagC)
{
  int tagn,tagCn,type,typeC;
  for (int n=0; n<npts; n++){
    tagn  = tag[n];
    tagCn = tagC[n];
    type  = bType[tagn];
    typeC = bType[tagCn];
    if (typeC > type) tag[n] = tagCn;
  }
}
