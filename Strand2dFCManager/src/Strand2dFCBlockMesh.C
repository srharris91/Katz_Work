/**
 * \brief
 * Implementation file for Class Strand2dFCBlockMesh.
 *
 * \author
 * Aaron Katz
 *
 * \version
 * 1.0
 *
 * \date
 * 11-22-2013
 */


#include "Strand2dFCBlockMesh.h"


// [Strand2dFCBlockMesh]
Strand2dFCBlockMesh::Strand2dFCBlockMesh()
{
  level = 0;
  meshOrder = 0;
  nSurfElem = 0;
  nSurfNode = 0;
  nBndNode = 0;
  nStrandNode = 0;
  nCompBd = 0;
}
// [Strand2dFCBlockMesh]


// [~Strand2dFCBlockMesh]
Strand2dFCBlockMesh::~Strand2dFCBlockMesh()
{
}
// [~Strand2dFCBlockMesh]


// [finalize]
void Strand2dFCBlockMesh::finalize()
{
  level = 0;
  meshOrder = 0;
  nSurfElem = 0;
  nSurfNode = 0;
  nBndNode = 0;
  nStrandNode = 0;
  nCompBd = 0;
  surfElem.deallocate();
  surfX.deallocate();
  bndNode.deallocate();
  surfElemTag.deallocate();
  bndNodeTag.deallocate();
  bndNodeNormal.deallocate();
  pointingVec.deallocate();
  clip.deallocate();
  strandX.deallocate();
}
// [finalize]


// get methods
const int& Strand2dFCBlockMesh::getMeshOrder(){return(meshOrder);}
const int& Strand2dFCBlockMesh::getNSurfElem(){return(nSurfElem);}
const int& Strand2dFCBlockMesh::getNSurfNode(){return(nSurfNode);}
const int& Strand2dFCBlockMesh::getNBndNode(){return(nBndNode);}
const int& Strand2dFCBlockMesh::getNStrandNode(){return(nStrandNode);}
const int& Strand2dFCBlockMesh::getNCompBd(){return(nCompBd);}
int* Strand2dFCBlockMesh::getSurfElem(){return(&surfElem(0,0));}
int* Strand2dFCBlockMesh::getBndNode(){return(&bndNode(0));}
double* Strand2dFCBlockMesh::getSurfX(){return(&surfX(0,0));}
double* Strand2dFCBlockMesh::getStrandX(){return(&strandX(0));}
int* Strand2dFCBlockMesh::getSurfElemTag(){return(&surfElemTag(0));}
int* Strand2dFCBlockMesh::getBndNodeTag(){return(&bndNodeTag(0));}
double* Strand2dFCBlockMesh::getBndNodeNormal(){return(&bndNodeNormal(0,0));}
double* Strand2dFCBlockMesh::getPointingVec(){return(&pointingVec(0,0));}
int* Strand2dFCBlockMesh::getClip(){return(&clip(0));}
