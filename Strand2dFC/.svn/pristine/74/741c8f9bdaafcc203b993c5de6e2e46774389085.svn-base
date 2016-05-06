#include "Strand2dFCSPTurbSA.h"
#include "Strand2dFCSPTurbSAInputRead.h"


void Strand2dFCSPTurbSA::inputRead(const string& inputFile)
{
  // hard coded values for this system layer
  nq          = 5; //rho,rho*u,rho*v,rho*e,rho*nu
  nqa         = 8; //P,u,v,T,nu,mu,k,dw
  ndim        = 2;
  ncomp       = 1;
  source      = 1; //add physical source terms

  // read values from input file
  int     nBpatchesT = 1000;
  int*    istate     = new int   [   ncomp];
  int*    itransport = new int   [   ncomp];
  int*    bTypeT     = new int   [   nBpatchesT];
  double* bValue     = new double[nq*nBpatchesT];
  
  strand2dfrspturbsainputread_(inputFile.size(),
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
  transport.getPrn(&Prn);
  transport.getPrnT(&PrnT);


  // initialize the Solution initialization class
  solution.initialize(isolution,
		      nq,
		      ndim,
		      inputFile,
		      &state,
		      &transport);


  // initialize the boundary condition class
  bType = new int[nBpatches];
  bc = new Strand2dFCSPTurbSABc*[nBpatches];
  for (int n=0; n<nBpatches; n++){
    bType[n] = bTypeT[n];
    if (bType[n] == 1) bc[n] = new Strand2dFCSPTurbSABcInviscidWall;
    if (bType[n] == 2) bc[n] = new Strand2dFCSPTurbSABcViscousWall;
    if (bType[n] == 3) bc[n] = new Strand2dFCSPTurbSABcFarField;
    if (bType[n] == 4) bc[n] = new Strand2dFCSPTurbSABcOutflow;
    bc[n]->initialize(nq,
		      nqa,
		      ndim,
		      ncomp,
		      viscous,
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
