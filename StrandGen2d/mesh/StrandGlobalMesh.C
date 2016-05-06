// File:        StrandGlobalMesh.C
// Package:     STRANDGEN
//              Parallel Infrastructure for Cartesian & Strand Solvers
// Release:     $Name:  $
// Revision:    $Revision:  $
// Modified:    $Date:  $
// Description: Data structure holding global mesh; partitions mesh and
//              creates StrandBlocks


#include "StrandGlobalMesh.h"


// Fortran interface definitions
extern "C" {
  void readinput_(const char*,
		  const int&,
		  char*,
		  const int&,
		  int&,
		  int&,
		  int&,
		  int&,
		  int&,
		  int&,
		  double&,
		  double&,
		  double&,
		  double&,
		  int&,
		  double&,
		  int&);
  void readmeshheader_(const char*,
		       const int&,
		       const int&,
		       int&,
		       int&,
		       int&,
		       int&,
		       int&,
		       int&,
		       int&,
		       int&);
  void readmeshdata_(const int&,
		     const int&,
		     const int&,
		     const int&,
		     const int&,
		     const int&,
		     const int&,
		     const int&,
		     const int&,
		     const int&,
		     int*,
		     int*,
		     int*,
		     int*,
		     int*,
		     int*,
		     int*,
		     double*,
		     double*,
		     double*,
		     double*);
  void strand1ddist_(const int&,
		     double*,
		     int&,
		     double&,
		     const int&,
		     double&,
		     double&);
  void findsharpcorners_(int&,
			 int&,
			 int&,
			 int&,
			 int*,
			 int*,
			 double*,
			 double&,
			 int&,
			 int&,
			 int*);
  void fillsharpcorners_(int&,
			 int&,
			 int&,
			 int&,
			 int&,
			 int*,
			 int*,
			 int*,
			 int*,
			 int*,
			 double*,
			 double*,
			 int&,
			 int*,
			 int&,
			 int&,
			 int*,
			 int*,
			 int*,
			 int*,
			 int*,
			 double*,
			 double*);
}


// static variable initialization
StrandGlobalMesh* StrandGlobalMesh::sgm_manager_instance = NULL;
bool    StrandGlobalMesh::d_initialized = false;
int     StrandGlobalMesh::d_nPartitions = 0;
int     StrandGlobalMesh::d_iwrite = 0;
int     StrandGlobalMesh::d_strandDist = 0;
int     StrandGlobalMesh::d_nPtsPerStrand = 0;
int     StrandGlobalMesh::d_nFringe = 0;
int     StrandGlobalMesh::d_trimType = 0;
int     StrandGlobalMesh::d_smax = 0;
double  StrandGlobalMesh::d_angle = 0.;
double  StrandGlobalMesh::d_stretchRatio = 0.;
double  StrandGlobalMesh::d_wallSpacing = 0.;
double  StrandGlobalMesh::d_strandLength = 0.;
double  StrandGlobalMesh::d_deltaSmooth = 0.;
int     StrandGlobalMesh::d_surfaceOnly = 0;
int     StrandGlobalMesh::d_nMesh = 0;
int     StrandGlobalMesh::d_nSurfNodes = 0;
int     StrandGlobalMesh::d_nSurfFaces = 0;
int     StrandGlobalMesh::d_nBndEdges = 0;
int     StrandGlobalMesh::d_nSurfPatches = 0;
int     StrandGlobalMesh::d_nEdgePatches = 0;
int     StrandGlobalMesh::d_nSharp = 0;
int*    StrandGlobalMesh::d_surfFaces = NULL;
int*    StrandGlobalMesh::d_bndEdges = NULL;
int*    StrandGlobalMesh::d_nodeClip = NULL;
int*    StrandGlobalMesh::d_faceClip = NULL;
int*    StrandGlobalMesh::d_faceTag = NULL;
int*    StrandGlobalMesh::d_bndTag = NULL;
int*    StrandGlobalMesh::d_sFlag = NULL;
double* StrandGlobalMesh::d_xSurf = NULL;
double* StrandGlobalMesh::d_pointingVec = NULL;
double* StrandGlobalMesh::d_xStrand = NULL;
double* StrandGlobalMesh::d_bndNormal = NULL;


// creation routine for singleton instance
void StrandGlobalMesh::createManager()
{
  if (!sgm_manager_instance) sgm_manager_instance = new StrandGlobalMesh();
}


// get the singleton instance
StrandGlobalMesh* StrandGlobalMesh::getManager()
{
  if (!sgm_manager_instance) createManager();
  return(sgm_manager_instance);
}


// close all open objects and de-allocate manager instance.
void StrandGlobalMesh::freeManager()
{
  if (sgm_manager_instance) delete sgm_manager_instance;
  sgm_manager_instance = ((StrandGlobalMesh*) NULL);
}


// protected Constructor
StrandGlobalMesh::StrandGlobalMesh()
{
  d_initialized = false;
  d_surfFaces   = NULL;
  d_bndEdges    = NULL;
  d_nodeClip    = NULL;
  d_faceClip    = NULL;
  d_faceTag     = NULL;
  d_bndTag      = NULL;
  d_sFlag       = NULL;
  d_xSurf       = NULL;
  d_pointingVec = NULL;
  d_xStrand     = NULL;
  d_bndNormal   = NULL;
}


// protected Destructor
StrandGlobalMesh::~StrandGlobalMesh()
{
  if (d_surfFaces  ) delete[] d_surfFaces;
  if (d_bndEdges   ) delete[] d_bndEdges;
  if (d_nodeClip   ) delete[] d_nodeClip;
  if (d_faceClip   ) delete[] d_faceClip;
  if (d_faceTag    ) delete[] d_faceTag;
  if (d_bndTag     ) delete[] d_bndTag;
  if (d_sFlag      ) delete[] d_sFlag;
  if (d_xSurf      ) delete[] d_xSurf;
  if (d_pointingVec) delete[] d_pointingVec;
  if (d_xStrand    ) delete[] d_xStrand;
  if (d_bndNormal  ) delete[] d_bndNormal;
  d_surfFaces   = NULL;
  d_bndEdges    = NULL;
  d_nodeClip    = NULL;
  d_faceClip    = NULL;
  d_faceTag     = NULL;
  d_bndTag      = NULL;
  d_sFlag       = NULL;
  d_xSurf       = NULL;
  d_pointingVec = NULL;
  d_xStrand     = NULL;
  d_bndNormal   = NULL;
}


// initialize the global mesh (read input and mesh files, allocate space for
// the global mesh and store it)
void StrandGlobalMesh::initialize(string& input_file){

  // read input file
  cout << "\nReading input from " << input_file << endl;
  int len = 80;
  int newlen;
  char* mfile = new char[len];
  readinput_(input_file.c_str(),
	     input_file.size(),
	     mfile,
	     len,
	     newlen,
	     d_nPartitions,
	     d_iwrite,
	     d_strandDist,
	     d_nPtsPerStrand,
	     d_nFringe,
	     d_stretchRatio,
	     d_wallSpacing,
	     d_strandLength,
	     d_deltaSmooth,
	     d_trimType,
	     d_angle,
	     d_smax);

  char* meshFile = new char[newlen];
  for (int n=0; n<newlen; n++) meshFile[n] = mfile[n];
  delete [] mfile;


  // read global mesh file header
  cout << "\nReading mesh from " << meshFile << endl;
  int iunit = 7;
  int nPtsPerStrand0;
  readmeshheader_(meshFile,
		  newlen,
		  iunit,
		  d_surfaceOnly,
		  d_nMesh,
		  d_nSurfNodes,
		  d_nSurfFaces,
		  d_nBndEdges,
		  d_nSurfPatches,
		  d_nEdgePatches,
		  nPtsPerStrand0);

  delete [] meshFile;
  if (d_surfaceOnly == 0) d_nPtsPerStrand = nPtsPerStrand0;


  // allocate space for the global mesh
  if (d_nSurfFaces > 0) d_surfFaces  = new int[2*d_nSurfFaces];
  if (d_nBndEdges  > 0){
    d_bndEdges  = new int   [  d_nBndEdges];
    d_bndTag    = new int   [  d_nBndEdges];
  }
  d_sFlag       = new int   [  d_nSurfNodes];
  d_nodeClip    = new int   [  d_nSurfNodes];
  d_faceClip    = new int   [  d_nSurfFaces];
  d_faceTag     = new int   [  d_nSurfFaces];
  d_xSurf       = new double[2*d_nSurfNodes];
  d_pointingVec = new double[2*d_nSurfNodes];
  d_xStrand     = new double[  d_nPtsPerStrand+1];
  if (d_nEdgePatches > 0) d_bndNormal = new double[2*d_nEdgePatches];

  // read mesh data
  readmeshdata_(iunit,
		d_surfaceOnly,
		d_nMesh,
		d_nSurfNodes,
		d_nSurfFaces,
		d_nBndEdges,
		d_nSurfPatches,
		d_nEdgePatches,
		d_nSharp,
		d_nPtsPerStrand,
		d_surfFaces,
		d_bndEdges,
		d_sFlag,
		d_nodeClip,
		d_faceClip,
		d_faceTag,
		d_bndTag,
		d_xSurf,
		d_pointingVec,
		d_xStrand,
		d_bndNormal);


  // locate sharp corners and add nodes and cells accordingly
  int  nNewNodes = 0;
  int* addNode = new int[d_nSurfNodes];
  for (int n; n<d_nSurfNodes; n++) addNode[n] = 1;

  if (d_angle > 0.){
    cout << "Detecting sharp corners..." << endl;

    findsharpcorners_(d_nSurfNodes,
		      d_nSurfFaces,
		      d_nSharp,
		      d_nBndEdges,
		      d_surfFaces,
		      d_bndEdges,
		      d_xSurf,
		      d_angle,
		      d_smax,
		      nNewNodes,
		      addNode);
  }

  int     nSurfNodesT  = d_nSurfNodes+nNewNodes;
  int     nSurfFacesT  = d_nSurfFaces+nNewNodes; //true in 2d only
  int*    surfFacesT   = new int   [2*nSurfFacesT];
  int*    nodeClipT    = new int   [  nSurfNodesT];
  int*    faceClipT    = new int   [  nSurfFacesT];
  int*    faceTagT     = new int   [  nSurfFacesT];
  int*    sFlagT       = new int   [  nSurfNodesT];
  double* xSurfT       = new double[2*nSurfNodesT];
  double* pointingVecT = new double[2*nSurfNodesT];

  fillsharpcorners_(d_nSurfNodes,
		    d_nSurfFaces,
		    d_nEdgePatches,
		    d_nSharp,
		    d_nBndEdges,
		    d_surfFaces,
		    d_bndEdges,
		    d_faceTag,
		    d_bndTag,
		    d_nodeClip,
		    d_xSurf,
		    d_bndNormal,
		    nNewNodes,
		    addNode,
		    nSurfNodesT,
		    nSurfFacesT,
		    surfFacesT,
		    nodeClipT,
		    faceClipT,
		    faceTagT,
		    sFlagT,
		    xSurfT,
		    pointingVecT);

  cout << nNewNodes << " nodes added" << endl;
  
  d_nSurfNodes = nSurfNodesT;
  d_nSurfFaces = nSurfFacesT;
  
  if (d_surfFaces  ) delete [] d_surfFaces;
  if (d_nodeClip   ) delete [] d_nodeClip;
  if (d_faceClip   ) delete [] d_faceClip;
  if (d_faceTag    ) delete [] d_faceTag;
  if (d_sFlag      ) delete [] d_sFlag;
  if (d_xSurf      ) delete [] d_xSurf;
  if (d_pointingVec) delete [] d_pointingVec;
  
  d_surfFaces   = new int   [2*d_nSurfFaces];
  d_nodeClip    = new int   [  d_nSurfNodes];
  d_faceClip    = new int   [  d_nSurfFaces];
  d_faceTag     = new int   [  d_nSurfFaces];
  d_sFlag       = new int   [  d_nSurfNodes];
  d_xSurf       = new double[2*d_nSurfNodes];
  d_pointingVec = new double[2*d_nSurfNodes];
  
  for (int n=0; n< 2*d_nSurfFaces; n++) d_surfFaces[n]   = surfFacesT[n];
  for (int n=0; n<   d_nSurfNodes; n++) d_nodeClip[n]    = nodeClipT[n];
  for (int n=0; n<   d_nSurfFaces; n++) d_faceClip[n]    = faceClipT[n];
  for (int n=0; n<   d_nSurfFaces; n++) d_faceTag[n]     = faceTagT[n];
  for (int n=0; n<   d_nSurfNodes; n++) d_sFlag[n]       = sFlagT[n];
  for (int n=0; n< 2*d_nSurfNodes; n++) d_xSurf[n]       = xSurfT[n];
  for (int n=0; n< 2*d_nSurfNodes; n++) d_pointingVec[n] = pointingVecT[n];
  
  if (surfFacesT  ) delete [] surfFacesT;
  if (nodeClipT   ) delete [] nodeClipT;
  if (faceClipT   ) delete [] faceClipT;
  if (faceTagT    ) delete [] faceTagT;
  if (sFlagT      ) delete [] sFlagT;
  if (xSurfT      ) delete [] xSurfT;
  if (pointingVecT) delete [] pointingVecT;
  

  // if mesh is only a surface mesh, find normal distribution base on
  // input file parameters
  if (d_surfaceOnly == 1){
    int nmax = 10000;
    double rmax;
    double* xs = new double[nmax];
    strand1ddist_(nmax,
		  xs,
		  d_nPtsPerStrand,
		  d_stretchRatio,
		  d_strandDist,
		  d_strandLength,
		  d_wallSpacing);
    delete [] d_xStrand;
    d_xStrand = new double[d_nPtsPerStrand+1];
    for (int n=0; n<d_nPtsPerStrand+1; n++) d_xStrand[n] = xs[n];
    delete [] xs;
  }


  // output strand distribution information
  cout.setf(ios::scientific);
  cout	<< "\nstrand distribution:";
  for (int n=0; n<d_nPtsPerStrand+1; n++)
    cout << "\n\t" << n << "\t" << *(d_xStrand+n);
  cout << "\n"
       << "\nwall spacing\t" << *(d_xStrand+1)
       << "\nend spacing\t"
       << d_xStrand[d_nPtsPerStrand]-d_xStrand[d_nPtsPerStrand-1]
       << "\nstretching ratio\t"
       << d_stretchRatio
       << "\n" << endl;


  // set initialization flag to true
  d_initialized = true;
}


const bool& StrandGlobalMesh::getInitialized()
{
  return(d_initialized);
}

const int& StrandGlobalMesh::getNPartitions()
{
   if (d_initialized) return(d_nPartitions);
   else{
     cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
     exit(0);
     return(d_nPartitions);
   }
}

const int& StrandGlobalMesh::getIwrite()
{
  if (d_initialized) return(d_iwrite);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_iwrite);
  }
}

const int& StrandGlobalMesh::getStrandDist()
{
  if (d_initialized) return(d_strandDist);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_strandDist);
  }
}

const int& StrandGlobalMesh::getNPtsPerStrand()
{
  if (d_initialized) return(d_nPtsPerStrand);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
     exit(0);
    return(d_nPtsPerStrand);
  }
}

const int& StrandGlobalMesh::getNFringe()
{
  if (d_initialized) return(d_nFringe);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_nFringe);
  }
}

const int& StrandGlobalMesh::getTrimType()
{
  if (d_initialized) return(d_trimType);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_trimType);
  }
}

const int& StrandGlobalMesh::getSmax()
{
  if (d_initialized) return(d_smax);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_smax);
  }
}

const double& StrandGlobalMesh::getAngle()
{
  if (d_initialized) return(d_angle);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_angle);
  }
}

const double& StrandGlobalMesh::getStretchRatio()
{
  if (d_initialized) return(d_stretchRatio);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_stretchRatio);
  }
}

const double& StrandGlobalMesh::getWallSpacing()
{
  if (d_initialized) return(d_wallSpacing);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_wallSpacing);
  }
}

const double& StrandGlobalMesh::getStrandLength()
{
  if (d_initialized) return(d_strandLength);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_strandLength);
  }
}

const double& StrandGlobalMesh::getDeltaSmooth()
{
  if (d_initialized) return(d_deltaSmooth);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_deltaSmooth);
  }
}

const int& StrandGlobalMesh::getSurfaceOnly()
{
  if (d_initialized) return(d_surfaceOnly);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_surfaceOnly);
  }
}

const int& StrandGlobalMesh::getNSurfNodes()
{
  if (d_initialized) return(d_nSurfNodes);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_nSurfNodes);
  }
}

const int& StrandGlobalMesh::getNSurfFaces()
{
  if (d_initialized) return(d_nSurfFaces);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_nSurfFaces);
  }
}

const int& StrandGlobalMesh::getNBndEdges()
{
  if (d_initialized) return(d_nBndEdges);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_nBndEdges);
  }
}

const int& StrandGlobalMesh::getNSurfPatches()
{
  if (d_initialized) return(d_nSurfPatches);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_nSurfPatches);
  }
}

const int& StrandGlobalMesh::getNEdgePatches()
{
  if (d_initialized) return(d_nEdgePatches);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_nEdgePatches);
  }
}

const int& StrandGlobalMesh::getNSharp()
{
  if (d_initialized) return(d_nSharp);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_nSharp);
  }
}

const int* StrandGlobalMesh::getSurfFaces()
{
  if (d_initialized) return(d_surfFaces);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_surfFaces);
  }
}

const int* StrandGlobalMesh::getBndEdges()
{
  if (d_initialized) return(d_bndEdges);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_bndEdges);
  }
}

const int* StrandGlobalMesh::getNodeClip()
{
  if (d_initialized) return(d_nodeClip);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_nodeClip);
  }
}

const int* StrandGlobalMesh::getFaceClip()
{
  if (d_initialized) return(d_faceClip);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_faceClip);
  }
}

const int* StrandGlobalMesh::getFaceTag()
{
  if (d_initialized) return(d_faceTag);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_faceTag);
  }
}

const int* StrandGlobalMesh::getBndTag()
{
  if (d_initialized) return(d_bndTag);
   else{
     cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
     exit(0);
     return(d_bndTag);
   }
}

const int* StrandGlobalMesh::getSFlag()
{
  if (d_initialized) return(d_sFlag);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_sFlag);
  }
}

const double* StrandGlobalMesh::getXSurf()
{
  if (d_initialized) return(d_xSurf);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_xSurf);
  }
}

const double* StrandGlobalMesh::getPointingVec()
{
  if (d_initialized) return(d_pointingVec);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_pointingVec);
  }
}

const double* StrandGlobalMesh::getXStrand()
{
  if (d_initialized) return(d_xStrand);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_xStrand);
  }
}

const double* StrandGlobalMesh::getBndNormal()
{
  if (d_initialized) return(d_bndNormal);
  else{
    cout << "StrandGlobalMesh: ERROR!! global mesh not initialized." << endl;
    exit(0);
    return(d_bndNormal);
  }
}

void StrandGlobalMesh::setNPartitions(int nPartitions)
{
  d_nPartitions = nPartitions;
}
