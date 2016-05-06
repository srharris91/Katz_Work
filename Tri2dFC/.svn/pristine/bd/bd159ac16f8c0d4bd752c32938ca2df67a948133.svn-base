#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::edgeExtractVis()
{
  int orderE=3-level; //order of local elements
  if      (orderE == 1){
    nee = 3;
    edgeE.allocate(nee,2);
    edgeE(0,0) = 0;
    edgeE(0,1) = 1;
    edgeE(1,0) = 0;
    edgeE(1,1) = 2;
    edgeE(2,0) = 1;
    edgeE(2,1) = 2;
  }
  else if (orderE == 2){
    nee = 9;
    edgeE.allocate(nee,2);
    edgeE(0,0) = 0;
    edgeE(0,1) = 3;
    edgeE(1,0) = 3;
    edgeE(1,1) = 1;
    edgeE(2,0) = 0;
    edgeE(2,1) = 5;
    edgeE(3,0) = 3;
    edgeE(3,1) = 5;
    edgeE(4,0) = 3;
    edgeE(4,1) = 4;
    edgeE(5,0) = 1;
    edgeE(5,1) = 4;
    edgeE(6,0) = 5;
    edgeE(6,1) = 4;
    edgeE(7,0) = 5;
    edgeE(7,1) = 2;
    edgeE(8,0) = 4;
    edgeE(8,1) = 2;
  }
  else if (orderE == 3){
    nee = 18;
    edgeE.allocate(nee,2);
    edgeE(0,0) = 0;
    edgeE(0,1) = 3;
    edgeE(1,0) = 3;
    edgeE(1,1) = 4;
    edgeE(2,0) = 4;
    edgeE(2,1) = 1;
    edgeE(3,0) = 0;
    edgeE(3,1) = 8;
    edgeE(4,0) = 3;
    edgeE(4,1) = 8;
    edgeE(5,0) = 3;
    edgeE(5,1) = 9;
    edgeE(6,0) = 4;
    edgeE(6,1) = 9;
    edgeE(7,0) = 4;
    edgeE(7,1) = 5;
    edgeE(8,0) = 1;
    edgeE(8,1) = 5;
    edgeE(9,0) = 8;
    edgeE(9,1) = 9;
    edgeE(10,0) = 9;
    edgeE(10,1) = 5;
    edgeE(11,0) = 8;
    edgeE(11,1) = 7;
    edgeE(12,0) = 9;
    edgeE(12,1) = 7;
    edgeE(13,0) = 9;
    edgeE(13,1) = 6;
    edgeE(14,0) = 5;
    edgeE(14,1) = 6;
    edgeE(15,0) = 7;
    edgeE(15,1) = 6;
    edgeE(16,0) = 7;
    edgeE(16,1) = 2;
    edgeE(17,0) = 6;
    edgeE(17,1) = 2;
  }
}
