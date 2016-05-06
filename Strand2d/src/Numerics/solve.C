#include "StrandBlockSolver.h"
#include "nbtrluC.h"
#include "nbtrbkC.h"


void StrandBlockSolver::solve(const int& step,
			      const int& pseudoStep,
			      const int& linearStep,
			      const int& mglevel,
			      const int& mode)
{
  // save the solution at the beginning of the pseudo step on the fine grid
  // for measuring RMS convergence
  if (mglevel == 0 && linearStep == 0){
    if (pseudoStep > 0){
      for (int k=0; k<nq; k++) rms(k) = 0.;
      for (int n=0; n<nFaces+nBedges; n++)
	for (int j=0; j<fClip(n)+1; j++)
	  for (int k=0; k<nq; k++)
	    rms(k) +=((q(k,j,n)-q0(k,j,n))*(q(k,j,n)-q0(k,j,n)));
      for (int k=0; k<nq; k++) rms(k) /=(rmsNorm(k)*rmsNorm(k));
    }
    for (int n=0; n<nFaces+nBedges; n++)
      for (int j=0; j<nPstr+2; j++)
	for (int k=0; k<nq; k++) q0(k,j,n) = q(k,j,n);
  }


  // initialize dq=0 and compute RHS at the beginning of the first linear step
  if (linearStep == 0){
    for (int n=0; n<nFaces+nBedges; n++)
      for (int j=0; j<nPstr+2; j++)
	for (int k=0; k<nq; k++) dq(k,j,n) = 0.;
    int sweep=0;
    computeRHS(step,
	       pseudoStep,
	       linearStep,
	       sweep,
	       mglevel,
	       mode);
  }


  // determine extents of forward or backward sweep
  int n1,n2,nn,m,i,l,jC,kk=1,jj,npts=nPstr+2;
  double a[nq];
  Array3D<double> de(nq,nq,npts);
  Array2D<double> rl(nq,npts);

  for (int sweep=0; sweep<2; sweep++){
    if      (sweep == 0){
      n1 = 0;
      n2 = nFaces+nBedges;
      nn = 1;
    }
    else{
      n1 = nFaces+nBedges-2;
      n2 = 0;
      nn =-1;
    }

    for (int n=n1; n!=n2; n+=nn){
      i  = gsMap(n);
      jC = fClip(i)+1;
      jj = fClip(i);
      if (standAlone == 1){
	jC = nPstr+2;
	jj = nPstr+1;
      }

      // move off-diagonal terms to RHS
      for (int j=0; j<jC; j++)
	for (int k=0; k<nq; k++) rl(k,j) = r(k,j,i);
      for (int m=ncsc(i); m<ncsc(i+1); m++){
	l = csc(m);
	for (int j=1; j<jC; j++){
	  matmul(nq,nq,kk,&bu(0,0,j,m),&dq(0,j,l),&a[0]);
	  for (int k=0; k<nq; k++) rl(k,j) += a[k];
	}}


      // invert tri-diagonal matrix for this strand stack and update
      for (int j=0; j<jC; j++)
	nbtrluC(0,jj,j,nq,&dm(0,0,j,i),&dd(0,0,j,i),&dp(0,0,j,i),
      		&rl(0,0),&de(0,0,0));
      nbtrbkC(0,jj,nq,&rl(0,0),&de(0,0,0));
      for (int j=0; j<jC; j++)
	for (int k=0; k<nq; k++) dq(k,j,i) = -rl(k,j);
    }}


  // update q and qa at the end of the last linear step
  if (linearStep == nLinearSteps-1){
    for (int n=0; n<nFaces+nBedges; n++)
      for (int j=0; j<nPstr+2; j++)
	for (int k=0; k<nq; k++) q(k,j,n) += dq(k,j,n);
    npts =(nFaces+nBedges)*(nPstr+2);
    sys->stepQAdd(npts,
		  &q(0,0,0),
		  &qa(0,0,0));
    nodalQFlag = 0;
    nodalQaFlag = 0;
    gradQFlag = 0;
    gradQaFlag = 0;
    limFlag = 0;
  }
  de.deallocate();
  rl.deallocate();
}
