#include "Strand2dFCManager.h"
#include "Strand2dFCManInputRead.h"
#include "strand1dDist.h"


void Strand2dFCManager::initialize(string& inputFile)
{
  nBlockMeshes  = 1; //for now just have one block
  nBlockSolvers = nBlockMeshes; //one solver per block


  // read input file
  int len=80;
  char meshFile[len];
  int meshOrders[len];
  int surfDist;
  int strandDist;
  int nStrandNodeG;
  double perturb;
  double wallSpacing;
  double strandLength;
  double stretchRatio;
  double deltaSmooth;
  strand2dfcmaninputread_(inputFile.size(),
			  inputFile.c_str(),
			  len,
			  meshFile,
			  nLevels,
			  &meshOrders[0],
			  surfDist,
			  strandDist,
			  perturb,
			  nStrandNodeG,
			  wallSpacing,
			  strandLength,
			  stretchRatio,
			  deltaSmooth,
			  iplotmesh);
  if (nLevels > meshOrders[0]){
    cout << "\n***Cannot set nLevels greater than finest meshOrder.***"
	 << endl;
  }


  // read mesh header
  cout << "\nReading global mesh from file " << meshFile << "\n" << endl;
  ifstream meshF;
  int nSurfElemG;
  int nSurfNodeG;
  int nBndNodeG;
  int nCompBdG;
  meshF.open(meshFile);
  meshF >> nSurfElemG >> nSurfNodeG >> nCompBdG >> nBndNodeG;


  // allocate data for the mesh, and read in from file
  int** surfElemG = new int*[nSurfElemG];
  Array1D<int> surfElemTagG(nSurfElemG);
  Array2D<double> surfXG(nSurfNodeG,2);
  Array1D<int> bndNodeG(nBndNodeG);
  Array1D<int> bndNodeTagG(nBndNodeG);
  Array2D<double> bndNodeNormalG(nBndNodeG,2);
  Array1D<double> strandXG(nStrandNodeG+1);
  int k,orderM;
  for (int n=0; n<nSurfElemG; n++){
    meshF >> orderM;
    surfElemG[n] = new int[orderM+2];
    surfElemG[n][0] = orderM;
    for (int m=1; m<orderM+2; m++) meshF >> surfElemG[n][m];
    meshF >> surfElemTagG(n);
  }
  for (int n=0; n<nSurfNodeG; n++) meshF >> surfXG(n,0) >> surfXG(n,1);
  for (int n=0; n<nBndNodeG; n++)
    meshF >> bndNodeG(n) >> bndNodeTagG(n)
	  >> bndNodeNormalG(n,0) >> bndNodeNormalG(n,1);
  meshF.close();


  // generate and potentially perturb 1-d strand distribution
  int nmax=10000,ntmp=nStrandNodeG;
  double rmax;
  double* xs = new double[nmax];
  strand1ddist_(nmax,
		xs,
		ntmp,
		stretchRatio,
		strandDist,
		strandLength,
		wallSpacing);
  nStrandNodeG++; //add end node
  for (int n=0; n<nStrandNodeG; n++) strandXG(n) = xs[n];
  Array1D<double> sx(nStrandNodeG);
  double dm,rn;
  for (int n=1; n<nStrandNodeG-1; n++){
    dm    = strandXG(n)-strandXG(n-1);
    dm    = max(dm,strandXG(n+1)-strandXG(n));
    rn    =((double)(rand()%100))/100.-.5;
    dm    = rn*perturb*dm;
    sx(n) = strandXG(n)+dm;
  }
  for (int n=1; n<nStrandNodeG-1; n++) strandXG(n) = sx(n);


  // allocate mesh blocks and initialize them
  blockMesh.allocate(nBlockMeshes,nLevels);
  for (int level=0; level<nLevels; level++)
    blockMesh(0,level).initialize(level,
				  meshOrders[level],
				  nSurfElemG,
				  nSurfNodeG,
				  nBndNodeG,
				  nStrandNodeG,
				  nCompBdG,
				  surfElemG,
				  bndNodeG,
				  surfXG,
				  strandXG,
				  surfElemTagG,
				  bndNodeTagG,
				  bndNodeNormalG);


  // clean up
  for (int n=0; n<nSurfElemG; n++){
    if (surfElemG[n]) delete [] surfElemG[n];
    surfElemG[n] = NULL;
  }
  delete [] surfElemG;
  surfElemG = NULL;
  delete [] xs;
  surfElemTagG.deallocate();
  surfXG.deallocate();
  bndNodeG.deallocate();
  bndNodeTagG.deallocate();
  bndNodeNormalG.deallocate();
  strandXG.deallocate();
  sx.deallocate();
}
