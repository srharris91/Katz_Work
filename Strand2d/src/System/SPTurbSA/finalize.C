#include "StrandSPTurbSA.h"


void StrandSPTurbSA::finalize()
{
  if (bc){
    for (int n=0; n<nBpatches; n++){
      if (bc[n]){
	bc[n]->finalize();
	delete bc[n];
	bc[n] = NULL;
      }
    }
    delete [] bc;
    bc = NULL;
  }

  state.finalize();
  transport.finalize();
  solution.finalize();
}
