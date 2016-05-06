#include "Tri2dFCManager.h"


void Tri2dFCManager::plot()
{
  if (iplotmesh != 0){
    ofstream meshF;
    meshF.open("mesh.vtk");
    meshF << "# vtk DataFile Version 3.0" << endl;
    meshF << "Global Mesh" << endl;
    meshF << "ASCII" << endl;
    meshF << "DATASET UNSTRUCTURED_GRID" << endl;
    meshF << "POINTS " << nNode << " float" << endl;
    for (int n=0; n<nNode; n++)
      meshF << x(n,0) << " " << x(n,1) << " " << 0. << endl;
    meshF << "CELLS " << nTri << " " << 4*nTri << endl;
    for (int n=0; n<nTri;  n++){
      meshF << 3 << " "
	    << tri[n][1] << " "
	    << tri[n][2] << " "
	    << tri[n][3] << endl;
    }
    meshF << "CELL_TYPES " << nTri << endl;
    for (int n=0; n<nTri;  n++) meshF << 5 << endl;
    meshF.close();

    cout << "\n *** mesh.vtk written in Paraview format. ***\n" << endl;
  }
}
