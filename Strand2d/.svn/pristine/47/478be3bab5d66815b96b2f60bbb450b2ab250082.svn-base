#include "StrandSPLam.h"
#include "StrandSPLamInputRead.h"


void StrandSPLam::inputRead(const string& inputFile)
{
  // hard coded values for this system layer
  nq          = 4; //rho,rho*u,rho*v,rho*e
  nqa         = 6; //P,u,v,T,mu,k
  ndim        = 2;
  ncomp       = 1;
  inviscid    = 1;
  dissipation = 1;
  source      = 0;

  // read values from input file
  int     nBpatchesT = 1000;
  int*    istate     = new int   [   ncomp     ];
  int*    itransport = new int   [   ncomp     ];
  int*    bType      = new int   [   nBpatchesT];
  double* bValue     = new double[nq*nBpatchesT];

  strandsplaminputread_(inputFile.size(),
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


  // initialize the Solution initialization class
  solution.initialize(isolution,
		      nq,
		      ndim,
		      inputFile,
		      &state,
		      &transport);


  // initialize the boundary condition class
  bc = new StrandSPLamBc*[nBpatches];
  for (int n=0; n<nBpatches; n++){
    if (bType[n] == 1) bc[n] = new StrandSPLamBcInviscidWall;
    if (bType[n] == 2) bc[n] = new StrandSPLamBcViscousWall;
    if (bType[n] == 3) bc[n] = new StrandSPLamBcInflow;
    if (bType[n] == 4) bc[n] = new StrandSPLamBcOutflow;
    if (bType[n] == 5) bc[n] = new StrandSPLamBcFarField;
    if (bType[n] == 6) bc[n] = new StrandSPLamBcDirichlet;
    if (bType[n] == 7) bc[n] = new StrandSPLamBcFrozen;
    bc[n]->initialize(nq,
		      nqa,
		      ndim,
		      ncomp,
		      gamma,
		      rGas,
		      bValue+n*nq,
		      &transport,
		      &solution);
  }

  delete [] istate;
  delete [] itransport;
  delete [] bType;
  delete [] bValue;
}
