#include "Tri2dFCBlockSolver.h"
#include "Tri2dFCSPLam.h"
#include "tri2dFCInputRead.h"


void Tri2dFCBlockSolver::input(const string& inputFile)
{
  int len=80;
  char systemType[len];
  int ordersT[len];
  tri2dfcinputread_(inputFile.size(),
		    inputFile.c_str(),
		    len,
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
		    nRKStages,
		    implicit,
		    nLevels,
		    &ordersT[0],
		    spacing,
		    gradMethod,
		    mgCycle,
		    limiter,
		    timeAcc,
		    dtUnsteady,
		    cfl,
		    vnn,
		    smooth,
		    convLimit,
		    relax,
		    systemType);

  orders.allocate(nLevels);
  for (int n=0; n<nLevels; n++) orders(n) = ordersT[n];


  // allocate system layer instance and read system layer inputs
  if (!strcmp(systemType,"SPLam"))
    sys = new Tri2dFCSPLam;
  else{
    cout << "\nSystem layer not recognized in input.C: Terminating." << endl;
    exit(0);
  }

  sys->inputRead(inputFile);


  // setup the system layer
  int    tmp = 500;
  int    iqgradT [tmp];
  int    iqagradT[tmp];
  double rmsNormT[tmp];
  double dlimT   [tmp];
  int    outputVarLengthT[tmp];
  string outputVarsT[500];
  sys->prepSetup(iPrint,
                 iTest,
                 iDebug,
                 tmp,
		 iSolnFile,
		 iResdFile,
		 iErrFile,
                 nq,
                 nqa,
                 ndim,
                 inviscid,
                 viscous,
                 source,
                 sourceMMS,
                 dissipation,
                 nCompBd,
                 &iqgradT[0],
                 &iqagradT[0],
                 &dlimT[0],
                 &rmsNormT[0],
		 nOutputVars,
		 &outputVarLengthT[0],
		 &outputVarsT[0]);

  nOutputScalars = 0;
  nOutputVectors = 0;
  outputVarLength.allocate(nOutputVars);
  outputVars.allocate(nOutputVars);
  for (int n=0; n<nOutputVars; n++){
    outputVarLength(n) = outputVarLengthT[n];
    outputVars(n) = outputVarsT[n];
    if      (outputVarLength(n) == 1) nOutputScalars++;
    else if (outputVarLength(n) == 3) nOutputVectors++;
    else{
      cout << "\n*** Incorrect outputVarLength in input.C ***" << endl;
      exit(0);
    }}

  nqGradQ = 0;
  for (int n=0; n<nq; n++) if (iqgradT[n] == 1) nqGradQ++;
  iqgrad.allocate(nqGradQ);
  nqGradQ = 0;
  for (int n=0; n<nq; n++)
    if (iqgradT[n] == 1) iqgrad(nqGradQ++) = n;

  nqaGradQa = 0;
  for (int n=0; n<nqa; n++) if (iqagradT[n] == 1) nqaGradQa++;
  iqagrad.allocate(nqaGradQa);
  nqaGradQa = 0;
  for (int n=0; n<nqa; n++)
    if (iqagradT[n] == 1) iqagrad(nqaGradQa++) = n;

  rmsNorm.allocate(nq);
  for (int n=0; n<nq; n++) rmsNorm(n) = rmsNormT[n];

  dlim.allocate(nq);
  for (int n=0; n<nq; n++) dlim(n) = dlimT[n];

  rka.allocate(nRKStages);
  rkb.allocate(nRKStages);
  if (nRKStages == 5){
    rka(0) = .25;
    rka(1) = 1./6.;
    rka(2) = .375;
    rka(3) = .5;
    rka(4) = 1.;
    rkb(0) = 1.;
    rkb(1) = 0.;
    rkb(2) = .56;
    rkb(3) = 0.;
    rkb(4) = .44;
  }
  else{
    cout << "\nCan only handle nRKStages = 5. Terminating" << endl;
    exit(0);
  }

  bdf.allocate(timeAcc+1);
  if      (timeAcc == 1){
    bdf(0) = 1.;
    bdf(1) =-1.;
  }
  else if (timeAcc == 2){
    bdf(0) = 1.5;
    bdf(1) =-2.;
    bdf(2) = .5;
  }
  else if (timeAcc == 3){
    bdf(0) = 11./6.;
    bdf(1) =-3.;
    bdf(2) = 1.5;
    bdf(3) =-1./3.;
  }
  else{
    cout << "\nCan only handle timeAcc = 1, 2, or 3. Terminating" << endl;
    exit(0);
  }
}
