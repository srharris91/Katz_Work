#include "Tri2dFCManager.h"


void Tri2dFCManager::convHistory(const int& step,
				 const int& pseudoStep,
				 bool& converged)
{
  cout.setf(ios::scientific);
  double* rms=t2dfcbs[0].getRms(); double rmsM=0.;
  for (int k=0; k<nq; k++){
    rms[k] = sqrt(rms[k]/(double)nDofs);
    if (rms[k] > rmsM) rmsM = rms[k];
  }
  if (rmsM < convLimit) converged = true;

  time = clock();
  double t=double(time-time0)/(double)CLOCKS_PER_SEC;
  cout << step << "\t" << pseudoStep << "\t" << t << "\t";
  for (int k=0; k<nq; k++) cout << rms[k] << "\t";
  cout << endl;

  if (cfile.is_open()){
    cfile << step << "\t" << pseudoStep << "\t" << t << "\t";
    for (int k=0; k<nq; k++) cfile << rms[k] << "\t";
    cfile << endl;
    if (converged) cfile.close();
  }

  int nPseudoStepsN=nPseudoSteps0;
  if (step > 0) nPseudoStepsN = nPseudoSteps;
  if (pseudoStep == nPseudoStepsN-1 || converged){
    double tS=double(time-timeS0)/(double)CLOCKS_PER_SEC;
    cout << "\nTime per pseudoStep for this physical step: "
	 << tS/(double)(pseudoStep+1) << "\n" << endl;
  }
}

