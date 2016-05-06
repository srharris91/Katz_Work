#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::sourceStencil()
{
  int ni,i1,i2,n1,n2;
  psp2S.allocate(nSurfNode+1);
  psp2S.set(0);

  if (surfOrder == 1 || surfOrder == 2) //mass lumped
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<meshOrder+1; i++) psp2S(surfElem(n,i)+1)++;

  /*
  else if (surfOrder == 2) //Galerkin
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<nElemEdge; i++){
	i1 = elemEdge(i,0);
	i2 = elemEdge(i,1);
	n1 = surfElem(n,i1);
	n2 = surfElem(n,i2);
	psp2S(n1+1) = psp2S(n1+1)+2;
	psp2S(n2+1) = psp2S(n2+1)+2;
      }
  */

  else if (surfOrder == 3) //corrected source by edge
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<meshOrder+1; i++)
	psp2S(surfElem(n,i)+1) = psp2S(surfElem(n,i)+1)+meshOrder+1;

  /*
  for (int n=0; n<nNode; n++)
    cout << n << " " << psp2S(n+1) << endl;
  exit(0);
  */

  for(int n=1; n<nSurfNode+1; n++) psp2S(n) += psp2S(n-1);
  npsp1S = psp2S(nSurfNode);
  psp1S.allocate(npsp1S,2);

  if (surfOrder == 1 || surfOrder == 2) //mass lumped
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<meshOrder+1; i++){
	ni                   = surfElem(n,i);
	psp1S(psp2S(ni)  ,0) = n;
	psp1S(psp2S(ni)++,1) = i;
      }

  /*
  if (surfOrder == 2) //Galerkin
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<nElemEdge; i++){
	i1 = elemEdge(i,0);
	i2 = elemEdge(i,1);
	n1 = surfElem(n,i1);
	n2 = surfElem(n,i2);
	psp1S(psp2S(n1)  ,0) = n;
	psp1S(psp2S(n1)++,1) = i1;
	psp1S(psp2S(n1)  ,0) = n;
	psp1S(psp2S(n1)++,1) = i2;
	psp1S(psp2S(n2)  ,0) = n;
	psp1S(psp2S(n2)++,1) = i1;
	psp1S(psp2S(n2)  ,0) = n;
	psp1S(psp2S(n2)++,1) = i2;
      }
  */

  else if (surfOrder == 3) //corrected source by edge
    for (int n=0; n<nSurfElem; n++)
      for (int i=0; i<meshOrder+1; i++){
	ni = surfElem(n,i);
	for (int m=0; m<meshOrder+1; m++){
	  psp1S(psp2S(ni)  ,0) = n;
	  psp1S(psp2S(ni)++,1) = m;
	}}

  for(int n=nSurfNode; n>0; n--) psp2S(n) = psp2S(n-1);
  psp2S(0) = 0;

  /*
  for (int n=0; n<nSurfNode; n++){
    cout << "\nSurfNode: " << n << " " << psp2S(n+1)-psp2S(n) << endl;
    for (int i=psp2S(n); i<psp2S(n+1); i++)
      cout << psp1S(i,0) << " " << psp1S(i,1) << endl;
  }
  exit(0);
  */
}
