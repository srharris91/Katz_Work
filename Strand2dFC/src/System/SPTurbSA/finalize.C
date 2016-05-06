#include "Strand2dFCSPTurbSA.h"


void Strand2dFCSPTurbSA::finalize()
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

  if (bType) delete [] bType;
  bType = NULL;

  state.finalize();
  transport.finalize();
  solution.finalize();
}
