#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::gradStencil()
{
  // create linked lists of points surrounding points
  int* esp2 = new int[nNode];
  for (int n=0; n<nNode; n++) esp2[n] = 0;
  for (int n=0; n<nElem; n++)
    for (int k=0; k<nne; k++) esp2[elem(n,k)]++;

  int** esp1 = new int*[nNode];
  for (int n=0; n<nNode; n++){
    if (esp2[n] > 0) esp1[n] = new int[esp2[n]];
    else esp1[n] = NULL;
  }

  int j;
  for (int n=0; n<nNode; n++) esp2[n] = 0;
  for (int n=0; n<nElem; n++)
    for (int k=0; k<nne; k++){
      j = esp2[elem(n,k)];
      esp1[elem(n,k)][j] = n;
      esp2[elem(n,k)]++;
    }

  /*
  for (int n=0; n<nNode; n++){
    cout << n << " ";
    for (int j=0; j<esp2[n]; j++) cout << esp1[n][j] << " ";
    cout << endl;
  }
  exit(0);
  */

  int i,m;
  Array1D<int> flag(nNode);
  psp2.allocate(nNode+1);
  psp2.set(0);
  flag.set(-1);
  for (int n=0; n<nNode; n++)
    for (int j=0; j<esp2[n]; j++){
      i = esp1[n][j];
      for (int k=0; k<nne; k++){
	m = elem(i,k);
	if (flag(m) != n){
	  psp2(n+1)++;
	  flag(m) = n;
	}}}

  /*
  for (int n=0; n<nNode; n++)
    cout << n << " " << psp2(n+1) << endl;
  exit(0);
  */

  for(int n=1; n<nNode+1; n++) psp2(n) += psp2(n-1);
  npsp1 = psp2(nNode);
  psp1.allocate(npsp1);

  flag.set(-1);
  for (int n=0; n<nNode; n++)
    for (int j=0; j<esp2[n]; j++){
      i = esp1[n][j];
      for (int k=0; k<nne; k++){
	m = elem(i,k);
	if (flag(m) != n){
	  psp1(psp2(n)++) = m;
	  flag(m)         = n;
	}}}

  for(int n=nNode; n>0; n--) psp2(n) = psp2(n-1);
  psp2(0) = 0;

  /*
  for (int n=0; n<nNode; n++){
    cout << "\nNode: " << n << " " << psp2(n+1)-psp2(n) << endl;
    for (int i=psp2(n); i<psp2(n+1); i++) cout << psp1(i) << endl;
  }
  exit(0);
  */


  // deallocate work arrays
  for (int n=0; n<nNode; n++)
    if (esp1[n]) delete [] esp1[n];
  delete [] esp1;
  delete [] esp2;
  flag.deallocate();
}
