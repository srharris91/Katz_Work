#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::edgeExtractVis()
{
  // flag all the nodes in the first element with their local surface
  // element and strand index numbers
  int n1;
  Array2D<int> flag(nNode,2);
  flag.set(-1);
  for (int i=0; i<meshOrder+1; i++){
    n1 = surfElem(0,i);
    for (int j=0; j<meshOrder+1; j++){
      flag(strandMap(n1,j),0) = i;
      flag(strandMap(n1,j),1) = j;
    }}


  // find all the edges for which both nodes have been flagged (these
  // are the edges used in the viscous discretization)
  int n2,n3,n4;
  nee = 0;
  for (int n=0; n<nEdge; n++){
    n1 = edge(n,0);
    n2 = edge(n,1);
    if (flag(n1,0) > -1 && flag(n2,0) > -1) nee++;
  }
  edgeE.allocate(nee,4,2);
  edgeE.set(-1);
  nee = 0;
  for (int n=0; n<nEdge-nEdgeBd; n++){
    n1 = edge(n,0);
    n2 = edge(n,1);
    n3 = edge(n,2);
    n4 = edge(n,3);
    if (flag(n1,0) > -1 && flag(n2,0) > -1){
      edgeE(nee,0,0) = flag(n1,0);
      edgeE(nee,0,1) = flag(n1,1);
      edgeE(nee,1,0) = flag(n2,0);
      edgeE(nee,1,1) = flag(n2,1);
      if (flag(n3,0) > -1){
	edgeE(nee,2,0) = flag(n3,0);
	edgeE(nee,2,1) = flag(n3,1);
      }
      if (flag(n4,0) > -1){
	edgeE(nee,3,0) = flag(n4,0);
	edgeE(nee,3,1) = flag(n4,1);
      }
      nee++;
    }}

  for (int n=nEdge-nEdgeBd; n<nEdge; n++){
    n1 = edge(n,0);
    n2 = edge(n,1);
    n3 = edge(n,2);
    if (flag(n1,0) > -1 && flag(n2,0) > -1){
      edgeE(nee,0,0) = flag(n1,0);
      edgeE(nee,0,1) = flag(n1,1);
      edgeE(nee,1,0) = flag(n2,0);
      edgeE(nee,1,1) = flag(n2,1);
      if (flag(n3,0) > -1){
	edgeE(nee,2,0) = flag(n3,0);
	edgeE(nee,2,1) = flag(n3,1);
      }
      nee++;
    }}

  /*
  for (int i=0; i<nee; i++)
    cout << i << "   "
	 << edgeE(i,0,0) << " "
	 << edgeE(i,0,1) << "   "
	 << edgeE(i,1,0) << " "
	 << edgeE(i,1,1) << "   "
	 << edgeE(i,2,0) << " "
	 << edgeE(i,2,1) << "   "  
	 << edgeE(i,3,0) << " "
	 << edgeE(i,3,1) << endl;
  */

  flag.deallocate();
}
