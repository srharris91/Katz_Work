#include "StrandBlockSolver.h"
#include "faceArea.h"
#include "cellCoord.h"
#include "volume.h"


void StrandBlockSolver::shiftTime(const int& step,
				  const int& mglevel)
{
  if (mglevel == 0){ // fine MG level
    cout << "\n*** Preparing to solve step " << step << " ***" << endl;
    cout << "physical time " << double(step)*dtUnsteady << endl;
    cout << "time step " << dtUnsteady << endl;


    // advance q, v, and x
    if (step == 1){
      for (int n=0; n<nFaces+nBedges; n++)
	for (int j=0; j<nPstr+2; j++){
	  for (int k=0; k<nq; k++){
	    q1(k,j,n) = q(k,j,n);
	    q2(k,j,n) = q(k,j,n);
	  }
	  v2(j,n) = v(j,n);
	  v1(j,n) = v(j,n);
	}
      for (int n=0; n<nNodes; n++)
	for (int j=0; j<nPstr+1; j++){
	  for (int k=0; k<ndim; k++){
	    x0(k,j,n) = x(k,j,n);
	    x1(k,j,n) = x(k,j,n);
	    x2(k,j,n) = x(k,j,n);
	  }}
    }
    else{
      for (int n=0; n<nFaces+nBedges; n++)
	for (int j=0; j<nPstr+2; j++){
	  for (int k=0; k<nq; k++){
	    q2(k,j,n) = q1(k,j,n);
	    q1(k,j,n) = q(k,j,n);
	  }
	  v2(j,n) = v1(j,n);
	  v1(j,n) = v(j,n);
	}
     for (int n=0; n<nNodes; n++)
	for (int j=0; j<nPstr+1; j++){
	  for (int k=0; k<ndim; k++){
	    x2(k,j,n) = x1(k,j,n);
	    x1(k,j,n) = x(k,j,n);
	    //x(k,j,n) = ?;
	  }}
    }


    // compute face areas, cell-centers, and volumes
    facearea_(pid,
	      ndim,
	      nFaces,
	      nGfaces,
	      nNodes,
	      nGnodes,
	      nEdges,
	      nBedges,
	      nPstr,
	      &face(0,0),
	      &edge(0,0),
	      &edgn(0),
	      &x(0,0,0),
	      &facs(0,0,0),
	      &facu(0,0,0));
    cellcoord_(pid,
	       ndim,
	       nFaces,
	       nGfaces,
	       nNodes,
	       nGnodes,
	       nEdges,
	       nBedges,
	       nPstr,
	       &face(0,0),
	       &edge(0,0),
	       &edgn(0),
	       &x(0,0,0),
	       &facs(0,0,0),
	       &facu(0,0,0),
	       &xc(0,0,0));
    volume_(pid,
	    ndim,
	    nFaces,
	    nGfaces,
	    nNodes,
	    nGnodes,
	    nEdges,
	    nBedges,
	    nPstr,
	    &face(0,0),
	    &edge(0,0),
	    &edgn(0),
	    &fClip(0),
	    &x(0,0,0),
	    &facs(0,0,0),
	    &facu(0,0,0),
	    &xc(0,0,0),
	    &v(0,0));


    // compute least squares coefficients
    if      (nodeVal == 1) lspVol();
    else if (nodeVal == 2) lspLS();
    else if (nodeVal == 3) lspMap();
    else{
      cout << "\nvalue of nodeVal not recognized" << endl;
      exit(0);
    }
  }


  else{ // coarse MG level
    // compute coarse level face areas and volumes
    coarseMetrics();
  }


  // face velocities
  cout << "\n*** face velocity computation not yet implemented ***" << endl;
}
