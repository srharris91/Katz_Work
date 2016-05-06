#include "StrandSPTurbSA.h"
#include "StrandSPTurbSAInputRead.h"


void StrandSPTurbSA::inputRead(const string& inputFile)
{
  // hard coded values for this system layer
  nq          = 5; //rho,rho*u,rho*v,rho*e,rho*nu_tilde
  nqa         = 8; //P,u,v,T,nu_tilde,mu_lam,k_lam,dw
  ndim        = 2;
  ncomp       = 1;
  inviscid    = 1;
  dissipation = 1;
  source      = 1;

  // read values from input file
  int     nBpatchesT = 1000;
  int*    istate     = new int   [   ncomp     ];
  int*    itransport = new int   [   ncomp     ];
  int*    bType      = new int   [   nBpatchesT];
  double* bValue     = new double[nq*nBpatchesT];

  strandspturbsainputread_(inputFile.size(),
			   inputFile.c_str(),
			   nq,
			   ncomp,
			   nBpatches,
			   nBpatchesT,
			   viscous,
			   sourceMMS,
			   isolution,
			   istate,
			   itransport,
			   bType,
			   bValue);


  // initialize the State class (Physics layer)
  state.initialize(ncomp,
		   istate,
		   inputFile);
  state.getGamma(&gamma);
  state.getRGas(&rGas);
  gm1  = gamma-1.;
  ggm1 = gamma/gm1;


  // initialize the Transport class (Physics layer)
  transport.initialize(ncomp,
		       itransport,
		       inputFile,
		       &state);
  transport.getPrnT(&prnT);


  // initialize the Solution initialization class
  solution.initialize(isolution,
		      nq,
		      ndim,
		      inputFile,
		      &state,
		      &transport);


  // initialize the boundary condition class
  bc = new StrandSPTurbSABc*[nBpatches];
  for (int n=0; n<nBpatches; n++){
    if (bType[n] == 1) bc[n] = new StrandSPTurbSABcInviscidWall;
    if (bType[n] == 2) bc[n] = new StrandSPTurbSABcViscousWall;
    if (bType[n] == 3) bc[n] = new StrandSPTurbSABcInflow;
    if (bType[n] == 4) bc[n] = new StrandSPTurbSABcOutflow;
    if (bType[n] == 5) bc[n] = new StrandSPTurbSABcFarField;
    if (bType[n] == 6) bc[n] = new StrandSPTurbSABcDirichlet;
    if (bType[n] == 7) bc[n] = new StrandSPTurbSABcFrozen;
    bc[n]->initialize(nq,
		      nqa,
		      ndim,
		      ncomp,
		      gamma,
		      rGas,
		      prnT,
		      bValue+n*nq,
		      &transport,
		      &solution);
  }

  delete [] istate;
  delete [] itransport;
  delete [] bType;
  delete [] bValue;
}
