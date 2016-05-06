#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::hessianStencil()
{
  // form prismatic elements from surface and strand elements
  int n1,j1,k,nElem=nSurfElem*nStrandElem;
  nne =(meshOrder+1)*(meshOrder+1);
  Array2D<int> elem(nElem,nne);
  nElem = 0;
  for (int n=0; n<nSurfElem; n++)
    for (int j=0; j<nStrandElem; j++){
      j1 = j*meshOrder;
      k  = 0;
      for (int i=0; i<meshOrder+1; i++){
	n1 = surfElem(n,i);
	for (int jj=0; jj<meshOrder+1; jj++)
	  elem(nElem,k++) = strandMap(n1,j1+jj);
      }
      nElem++;
    }

  /*
  for (int n=0; n<nElem; n++){
    cout << n << "  ";
    for (int i=0; i<nne; i++) cout << elem(n,i) << " ";
    cout << endl;
  }
  exit(0);
  */


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
  pspH2.allocate(nNode+1);
  pspH2.set(0);
  flag.set(-1);
  for (int n=0; n<nNode; n++)
    for (int j=0; j<esp2[n]; j++){
      i = esp1[n][j];
      for (int k=0; k<nne; k++){
	m = elem(i,k);
	if (flag(m) != n){
	  pspH2(n+1)++;
	  flag(m) = n;
	}}}

  /*
  for (int n=0; n<nNode; n++)
    cout << n << " " << pspH2(n+1) << endl;
  exit(0);
  */

  for(int n=1; n<nNode+1; n++) pspH2(n) += pspH2(n-1);
  npspH1 = pspH2(nNode);
  pspH1.allocate(npspH1);

  flag.set(-1);
  for (int n=0; n<nNode; n++)
    for (int j=0; j<esp2[n]; j++){
      i = esp1[n][j];
      for (int k=0; k<nne; k++){
	m = elem(i,k);
	if (flag(m) != n){
	  pspH1(pspH2(n)++) = m;
	  flag(m)           = n;
	}}}

  for(int n=nNode; n>0; n--) pspH2(n) = pspH2(n-1);
  pspH2(0) = 0;

  /*
  for (int n=0; n<nNode; n++){
    cout << "\nNode: " << n << " " << pspH2(n+1)-pspH2(n) << endl;
    for (int i=pspH2(n); i<pspH2(n+1); i++) cout << pspH1(i) << endl;
  }
  exit(0);
  */


  // deallocate work arrays
  for (int n=0; n<nNode; n++)
    if (esp1[n]) delete [] esp1[n];
  delete [] esp1;
  delete [] esp2;
  flag.deallocate();
  elem.deallocate();
}
