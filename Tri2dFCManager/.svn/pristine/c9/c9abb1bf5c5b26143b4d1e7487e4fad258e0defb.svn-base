#include "Tri2dFCManager.h"
#include "tri2dFCManInputRead.h"


void Tri2dFCManager::initialize(string& inputFile)
{
  // read input file
  int len=80;
  char meshFile[len];
  tri2dfcmaninputread_(inputFile.size(),
		       inputFile.c_str(),
		       iplotmesh,
		       len,
		       meshFile);


  // read mesh header
  cout << "\nReading global mesh from file " << meshFile << "\n" << endl;
  ifstream meshF;
  meshF.open(meshFile);
  meshF >> nTri >> nNode >> nCompBd >> nEdgeBd;
  cout << "\nMesh statistics: " << endl
       << "Number of triangles: " << nTri << endl
       << "Number of Nodes: " << nNode << endl
       << "Number of boundary comp.: " << nCompBd << endl
       << "Number of boundary edges: " << nEdgeBd << "\n" << endl;


  // allocate mesh data
  tri = new int*[nTri];
  x.allocate(nNode,2);
  edgeBd.allocate(nEdgeBd,3);


  // read mesh data
  int k,ord;
  for (int n=0; n<nTri;  n++){
    meshF >> ord;
    k         =(ord+2)*(ord+1)/2;
    tri[n]    = new int[k+1];
    tri[n][0] = ord;
    for (int m=1; m<k+1; m++) meshF >> tri[n][m];
  }
  for (int n=0; n<nNode; n++) meshF >> x(n,0) >> x(n,1);
  for (int n=0; n<nEdgeBd; n++)
    meshF >> edgeBd(n,0) >> edgeBd(n,1) >> edgeBd(n,2);

  meshF.close();
}
