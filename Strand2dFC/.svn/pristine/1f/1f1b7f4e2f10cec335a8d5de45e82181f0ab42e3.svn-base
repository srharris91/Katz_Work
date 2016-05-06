#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::gradStencil()
{
  // create linked lists of surface elements surrounding surface points
  int* surfEsp2 = new int[nSurfNode];
  for (int n=0; n<nSurfNode; n++) surfEsp2[n] = 0;
  for (int n=0; n<nSurfElem; n++)
    for (int k=0; k<meshOrder+1; k++) surfEsp2[surfElem(n,k)]++;

  int** surfEsp1 = new int*[nSurfNode];
  for (int n=0; n<nSurfNode; n++){
    if (surfEsp2[n] > 0) surfEsp1[n] = new int[surfEsp2[n]];
    else surfEsp1[n] = NULL;
  }

  int j;
  for (int n=0; n<nSurfNode; n++) surfEsp2[n] = 0;
  for (int n=0; n<nSurfElem; n++)
    for (int k=0; k<meshOrder+1; k++){
      j = surfEsp2[surfElem(n,k)];
      surfEsp1[surfElem(n,k)][j] = n;
      surfEsp2[surfElem(n,k)]++;
    }

  /*
  for (int n=0; n<nSurfNode; n++){
    cout << n << " ";
    for (int j=0; j<surfEsp2[n]; j++) cout << surfEsp1[n][j] << " ";
    cout << endl;
  }
  exit(0);
  */


  int i,m,n1,n2;
  Array1D<int> flag(nSurfNode);
  psp2.allocate(nSurfNode+1);
  psp2.set(0);
  flag.set(-1);
  for (int n=0; n<nSurfNode; n++)
    for (int l=0; l<surfEsp2[n]; l++){
      i = surfEsp1[n][l];
      for (int k=0; k<meshOrder+1; k++){
	m = surfElem(i,k);
	if (flag(m) != n){
	  psp2(n+1)++;
	  flag(m) = n;
	}}}

  /*
  for (int n=0; n<nNode; n++)
    cout << n << " " << psp2(n+1) << endl;
  exit(0);
  */

  for(int n=1; n<nSurfNode+1; n++) psp2(n) += psp2(n-1);
  npsp1 = psp2(nSurfNode);
  psp1.allocate(npsp1);

  flag.set(-1);
  for (int n=0; n<nSurfNode; n++)
    for (int l=0; l<surfEsp2[n]; l++){
      i = surfEsp1[n][l];
      for (int k=0; k<meshOrder+1; k++){
	m = surfElem(i,k);
	if (flag(m) != n){
	  psp1(psp2(n)++) = m;
	  flag(m)         = n;
	}}}

  for(int n=nSurfNode; n>0; n--) psp2(n) = psp2(n-1);
  psp2(0) = 0;

  /*
  for (int n=0; n<nSurfNode; n++){
    cout << "\nSurfNode: " << n << " " << psp2(n+1)-psp2(n) << endl;
    for (int i=psp2(n); i<psp2(n+1); i++) cout << psp1(i) << endl;
  }
  exit(0);
  */


  // deallocate work arrays
  for (int n=0; n<nSurfNode; n++)
    if (surfEsp1[n]) delete [] surfEsp1[n];
  delete [] surfEsp1;
  delete [] surfEsp2;
  flag.deallocate();
}
