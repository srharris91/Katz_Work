#include "StrandBlockSolver.h"


void StrandBlockSolver::order()
{
  // initialize gsMap
  for (int n=1; n<nFaces+nBedges; n++) gsMap(n) = -1;


  // find starting cell
  int cc=0,j=nEdges-1;
  if (nBedges > 0) cc = edge(1,j);
  gsMap(0) = cc;
  int flag[nFaces+nBedges];
  for (int n=0; n<nFaces+nBedges; n++) flag[n] = 0;
  flag[cc] = 1;


  // loop through cells, finding adjacent cells until all cells are gone
  for (int n=1; n<nFaces+nBedges; n++){
    for (int m=ncsc(cc); m<ncsc(cc+1); m++){
      j = csc(m);
      if (flag[j] == 0){
	cc       = j;
	flag[cc] = 1;
	gsMap(n) = cc;
	break;
      }}}


  for (int n=0; n<nFaces+nBedges; n++)
    if (gsMap(n) < 0){
      cout << "\ngsMap procedure failure. Terminating." << endl;
      exit(0);
    }
}
