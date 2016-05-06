#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::reorderNodes()
{
  // resolve boundary types at corners with the system layer
  int k;
  Array2D<int> nflag(nNode,2);
  nflag.set(-1);
  for (int n=nEdge-nEdgeBd; n<nEdge; n++)
    for (int j=0; j<2; j++){
      k  = edge(n,j);
      if (nflag(k,0) >= 0) sys->bcConflict(1,&nflag(k,0),&edge(n,3));
      else nflag(k,0) = edge(n,3);
      //nflag(k,0) = edge(n,3);
    }


  // form boundary node array
  k = 0;
  for (int n=0; n<nNode; n++)
    if (nflag(n,0) < 0) nflag(n,1) = k++;
  nNodeBd = 0;
  for (int n=0; n<nNode; n++)
    if (nflag(n,0) >= 0){
      nflag(n,1) = k++;
      nNodeBd++;
    }

  if (k != nNode){
    cout << "\n*** problem re-ordering nodes in connectivity.C ***" << endl;
    exit(0);
  }

  /*
  for (int n=0; n<nNode; n++)
    cout << n << " " << nflag(n,1) << endl;
  exit(0);
  */

  nodeBd.allocate(nNodeBd);
  int m=0;
  for (int n=0; n<nNode; n++)
    if (nflag(n,0) >= 0) nodeBd(m++) = nflag(n,0);

  /*
  m = 0;
  for (int n=nNode-nNodeBd; n<nNode; n++){
    cout << n << " " << m << " " << nodeBd(m) << endl;
    m++;
  }
  */

  for (int n=0; n<nElem; n++)
    for (int j=0; j<nne; j++) elem(n,j) = nflag(elem(n,j),1);

  for (int n=0; n<nTri; n++)
    for (int j=0; j<3; j++) tri(n,j) = nflag(tri(n,j),1);

  Array2D<double> xT(nNode,2);
  for (int n=0; n<nNode; n++)
    for (int j=0; j<2; j++) xT(nflag(n,1),j) = x(n,j);
  for (int n=0; n<nNode; n++)
    for (int j=0; j<2; j++) x(n,j) = xT(n,j);

  for (int n=0; n<nEdgeBd; n++)
    for (int j=0; j<2; j++) edgeBd(n,j) = nflag(edgeBd(n,j),1);

  for (int n=0; n<nEdge-nEdgeBd; n++)
    for (int j=0; j<4; j++) edge(n,j) = nflag(edge(n,j),1);

  for (int n=nEdge-nEdgeBd; n<nEdge; n++)
    for (int j=0; j<3; j++) edge(n,j) = nflag(edge(n,j),1);


  // deallocate work arrays
  nflag.deallocate();
  xT.deallocate();

  /*
  for (int n=0; n<nEdge; n++)
    cout << n << " "
	 << edge(n,0) << " " << edge(n,1) << " "
	 << edge(n,2) << " " << edge(n,3) << " "
	 << edge(n,4) << " " << edge(n,5) << endl;
  */
}
