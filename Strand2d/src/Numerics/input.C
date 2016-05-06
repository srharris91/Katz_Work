#include "StrandBlockSolver.h"
#include "StrandSPLam.h"
#include "StrandSPTurbSA.h"
#include "strandInputRead.h"


void StrandBlockSolver::input(const string& inputFile)
{
  int len=80;
  char systemType[len];
  strandinputread_(inputFile.size(),
		   inputFile.c_str(),
		   len,
		   systemType,
		   iPrint,
		   iTest,
		   iDebug,
		   iConvFile,
		   iSolnFile,
		   iResdFile,
		   iErrFile,
		   iSurfFile,
		   standAlone,
		   restartStep,
		   nRestart,
		   nOutput,
		   nSteps,
		   nPseudoSteps,
		   nPseudoSteps0,
		   nLinearSteps,
		   nLevels,
		   mgCycle,
		   gradient,
		   nodeVal,
		   perturb,
		   limiter,
		   dtUnsteady,
		   cflLinear,
		   cfl0,
		   vnn0,
		   nRamp,
		   brelax,
		   gradClip,
		   cfl,
		   vnn,
		   convLimit,
		   relax,
		   coarseDis);

  if (!strcmp(systemType,"SPLam"))
    sys = new StrandSPLam;
  else if (!strcmp(systemType,"SPTurbSA"))
    sys = new StrandSPTurbSA;
  else{
    cout << "\nSystem layer not recognized in input.C: Terminating." << endl;
    exit(0);
  }

  sys->inputRead(inputFile);

  int tmp = 500;
  Array1D<int>    iqgradT (tmp);
  Array1D<int>    iqagradT(tmp);
  Array1D<double> rmsNormT(tmp);
  Array1D<double> dlimT   (tmp);

  sys->prepSetup(iPrint,
		 iTest,
		 iDebug,
		 tmp,
		 nq,
		 nqa,
		 ndim,
		 inviscid,
		 viscous,
		 source,
		 sourceMMS,
		 dissipation,
		 nBpatches,
		 &iqgradT(0),
		 &iqagradT(0),
		 &dlimT(0),
		 &rmsNormT(0));

  nqGradQ = 0;
  for (int n=0; n<nq; n++) if (iqgradT(n) == 1) nqGradQ++;
  iqgrad.allocate(nqGradQ);
  nqGradQ = 0;
  for (int n=0; n<nq; n++)
    if (iqgradT(n) == 1){
      iqgrad(nqGradQ) = n;
      nqGradQ++;
    }

  nqaGradQa = 0;
  for (int n=0; n<nqa; n++) if (iqagradT(n) == 1) nqaGradQa++;
  iqagrad.allocate(nqaGradQa);
  nqaGradQa = 0;
  for (int n=0; n<nqa; n++)
    if (iqagradT(n) == 1){
      iqagrad(nqaGradQa) = n;
      nqaGradQa++;
    }

  rmsNorm.allocate(nq);
  for (int n=0; n<nq; n++) rmsNorm(n) = rmsNormT(n);

  dlim.allocate(nq);
  for (int n=0; n<nq; n++) dlim(n) = dlimT(n);

  iqgradT.deallocate();
  iqagradT.deallocate();
  rmsNormT.deallocate();
  dlimT.deallocate();
}
