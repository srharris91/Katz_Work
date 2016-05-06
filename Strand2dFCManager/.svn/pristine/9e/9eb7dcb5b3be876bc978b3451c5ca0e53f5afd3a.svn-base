#include "Strand2dFCManager.h"


void Strand2dFCBlockMesh::plot()
{
  // create a fully unstructured grid from the strand grid
  // (plot only the finest level at this point)
  int nNode=nSurfNode*nStrandNode;
  int nQuad=nSurfElem*meshOrder*(nStrandNode-1);
  Array2D<double> x(nNode,2);
  Array2D<int> map(nSurfNode,nStrandNode);
  Array2D<int> quad(nQuad,4);
  Array2D<int> surfCell(meshOrder,2);
  if (meshOrder == 1){
    surfCell(0,0) = 0;
    surfCell(0,1) = 1;
  }
  else if (meshOrder == 2){
    surfCell(0,0) = 0;
    surfCell(0,1) = 2;
    surfCell(1,0) = 2;
    surfCell(1,1) = 1;
  }
  else if (meshOrder == 3){
    surfCell(0,0) = 0;
    surfCell(0,1) = 2;
    surfCell(1,0) = 2;
    surfCell(1,1) = 3;
    surfCell(2,0) = 3;
    surfCell(2,1) = 1;
  }
  else if (meshOrder == 4){
    surfCell(0,0) = 0;
    surfCell(0,1) = 2;
    surfCell(1,0) = 2;
    surfCell(1,1) = 3;
    surfCell(2,0) = 3;
    surfCell(2,1) = 4;
    surfCell(3,0) = 4;
    surfCell(3,1) = 1;
  }
  
  double x0,y0,nx,ny;
  nNode = 0;
  for (int n=0; n<nSurfNode; n++){
    x0 = surfX(n,0);
    y0 = surfX(n,1);
    nx = pointingVec(n,0);
    ny = pointingVec(n,1);
    for (int j=0; j<nStrandNode; j++){
      x(nNode,0) = x0+nx*strandX(j);
      x(nNode,1) = y0+ny*strandX(j);
      map(n,j)   = nNode++;
    }}
  
  int ord,n1,n2;
  nQuad = 0;
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<meshOrder; i++){
      n1 = surfElem(n,surfCell(i,0));
      n2 = surfElem(n,surfCell(i,1));
      for (int j=0; j<nStrandNode-1; j++){
	quad(nQuad  ,0) = map(n1,j  );
	quad(nQuad  ,1) = map(n2,j  );
	quad(nQuad  ,2) = map(n2,j+1);
	quad(nQuad++,3) = map(n1,j+1);
      }}
  
  
  // plot the fully unstructured grid
  stringstream a;
  a.str("");
  a.clear();
  a << "meshLevel" << level << ".vtk";
  ofstream meshF;
  meshF.open (a.str().c_str());
  meshF << "# vtk DataFile Version 3.0" << endl;
  meshF << "Global Mesh" << endl;
  meshF << "ASCII" << endl;
  meshF << "DATASET UNSTRUCTURED_GRID" << endl;
  meshF << "POINTS " << nNode << " float" << endl;
  for (int n=0; n<nNode; n++)
    meshF << x(n,0) << " " << x(n,1) << " " << 0. << endl;
  meshF << "CELLS " << nQuad << " " << 5*nQuad << endl;
  for (int n=0; n<nQuad; n++)
    meshF << 4 << " "
	  << quad(n,0) << " "
	  << quad(n,1) << " "
	  << quad(n,2) << " "
	  << quad(n,3) << endl;
  meshF << "CELL_TYPES " << nQuad << endl;
  for (int n=0; n<nQuad; n++) meshF << 9 << endl;
  meshF.close();
  
  x.deallocate();
  map.deallocate();
  quad.deallocate();
  surfCell.deallocate();
  
  cout << "\n *** " << a.str() << " written in Paraview format. ***\n" << endl;
}
