#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::edgeExtract()
{
  // form local edges in the standard element
  nElemEdge = meshOrder;
  elemEdge.allocate(nElemEdge,2);
  if (meshOrder == 1){
    elemEdge(0,0) = 0;
    elemEdge(0,1) = 1;
  }
  else if (meshOrder == 2){
    elemEdge(0,0) = 0;
    elemEdge(0,1) = 2;
    elemEdge(1,0) = 2;
    elemEdge(1,1) = 1;
  }
  else if (meshOrder == 3){
    elemEdge(0,0) = 0;
    elemEdge(0,1) = 2;
    elemEdge(1,0) = 2;
    elemEdge(1,1) = 3;
    elemEdge(2,0) = 3;
    elemEdge(2,1) = 1;
  }
  else if (meshOrder == 4){
    elemEdge(0,0) = 0;
    elemEdge(0,1) = 2;
    elemEdge(1,0) = 2;
    elemEdge(1,1) = 3;
    elemEdge(2,0) = 3;
    elemEdge(2,1) = 4;
    elemEdge(3,0) = 4;
    elemEdge(3,1) = 1;
  }

  /*
  for (int n=0; n<nElemEdge; n++)
    cout << n << " " << elemEdge(n,0) << " " << elemEdge(n,1) << endl;
  exit(0);
  */

  // form global list of edges
  // edge(n,0) = left node on edge n
  // edge(n,1) = right node on edge n
  // edge(n,2) = left limiter extension node on edge n
  // edge(n,3) = right limiter extension node on edge n
  nSurfEdge = nSurfElem*meshOrder;
  surfEdge.allocate(nSurfEdge,4);//NOTE: THIS SHOULD BE IN LOCAL ELEMENTS
  nSurfEdge = 0;
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<meshOrder; i++){
      surfEdge(nSurfEdge  ,0) = surfElem(n,elemEdge(i,0));
      surfEdge(nSurfEdge++,1) = surfElem(n,elemEdge(i,1));
    }


  // form edge extension nodes, setting the extension node to the opposite
  // edge node for edges that impinge on the boundary
  int n1,n2;
  Array2D<int> psp(nSurfNode,2);
  psp.set(-1);
  for (int n=0; n<nSurfEdge; n++){
    n1        = surfEdge(n,0);
    n2        = surfEdge(n,1);
    psp(n1,1) = n2;
    psp(n2,0) = n1;
  }
  for (int n=0; n<nSurfEdge; n++){
    n1            = surfEdge(n,0);
    n2            = surfEdge(n,1);
    if (psp(n1,0) == -1) surfEdge(n,2) = n2;
    else surfEdge(n,2) = psp(n1,0);
    if (psp(n2,1) == -1) surfEdge(n,3) = n1;
    else surfEdge(n,3) = psp(n2,1);
  }
  /*
  for (int n=0; n<nSurfEdge; n++)
    cout << n << " "
	 << surfEdge(n,2) << " "
	 << surfEdge(n,0) << " "
	 << surfEdge(n,1) << " "
	 << surfEdge(n,3) << " "
	 << endl;
  */


  // determine if boundary nodes are on the right or left side of the domain
  bndSign.allocate(nBndNode);
  bndElem.allocate(nBndNode,2);
  Array1D<int> flag(nSurfNode);
  flag.set(-1);
  for (int n=0; n<nBndNode; n++) flag(bndNode(n)) = n;
  for (int n=0; n<nSurfElem; n++){
    n1 = flag(surfElem(n,0));
    n2 = flag(surfElem(n,1));
    if (n1 > -1) bndSign(n1) = -1.; //left boundary
    if (n2 > -1) bndSign(n2) =  1.; //right boundary
    for (int i=0; i<meshOrder+1; i++){
      n1 = flag(surfElem(n,i));
      if (n1 > -1){
	bndElem(n1,0) = n;
	bndElem(n1,1) = i;
      }
    }}
  //for (int n=0; n<nBndNode; n++)
  //cout << bndNode(n) << " " << bndSign(n) << endl;


  // clean up
  psp.deallocate();
  flag.deallocate();
}
