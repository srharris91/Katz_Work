#include "Tri2dFCSPLam.h"
#include "Tri2dFCSPLamInputRead.h"


void Tri2dFCSPLam::inputRead(const string& inputFile)
{
  // hard coded values for this system layer
  nq          = 4; //rho,rho*u,rho*v,rho*e
  nqa         = 6; //P,u,v,T,mu,k
  ndim        = 2;
  ncomp       = 1;
  source      = 0;

  // read values from input file
  int     nBpatchesT = 1000;
  int*    istate     = new int   [   ncomp];
  int*    itransport = new int   [   ncomp];
  int*    bTypeT     = new int   [   nBpatchesT];
  double* bValue     = new double[nq*nBpatchesT];
  
  tri2dfrsplaminputread_(inputFile.size(),
			 inputFile.c_str(),
			 nq,
			 ncomp,
			 nBpatches,
			 nBpatchesT,
			 inviscid,
			 dissipation,
			 viscous,
			 sourceMMS,
			 isolution,
			 istate,
			 itransport,
			 bTypeT,
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
  nBpatches++; //add a "nothing" condition in the case of sharp corners
  bTypeT[nBpatches-1] = 10;
  bType = new int[nBpatches];
  for (int n=0; n<nq; n++) bValue[(nBpatches-1)*nq+n] = 0.;
  bc = new Tri2dFCSPLamBc*[nBpatches];
  for (int n=0; n<nBpatches; n++){
    bType[n] = bTypeT[n];
    if (bType[n] == 1 ) bc[n] = new Tri2dFCSPLamBcInviscidWall;
    if (bType[n] == 5 ) bc[n] = new Tri2dFCSPLamBcViscousWall;
    if (bType[n] == 3 ) bc[n] = new Tri2dFCSPLamBcInflow;
    if (bType[n] == 2 ) bc[n] = new Tri2dFCSPLamBcOutflow;
    if (bType[n] == 4 ) bc[n] = new Tri2dFCSPLamBcFarField;
    if (bType[n] == 6 ) bc[n] = new Tri2dFCSPLamBcDirichlet;
    if (bType[n] == 7 ) bc[n] = new Tri2dFCSPLamBcFrozen;
    if (bType[n] == 8 ) bc[n] = new Tri2dFCSPLamBcSymmetry;
    if (bType[n] == 9 ) bc[n] = new Tri2dFCSPLamBcPeriodic;
    if (bType[n] == 10) bc[n] = new Tri2dFCSPLamBcNothing;
    bc[n]->initialize(nq,
		      nqa,
		      ndim,
		      ncomp,
		      gamma,
		      rGas,
		      bValue+n*nq,
		      &state,
		      &transport,
		      &solution);
  }

  delete [] istate;
  delete [] itransport;
  delete [] bTypeT;
  delete [] bValue;
}
