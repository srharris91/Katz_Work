#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::faceArea()
{
  int p1,p2,n1,n2;
  double xa,ya,xb,yb;
  for (int n=0; n<nEdge-nEdgeBd; n++){ //interior face areas
    n1        = edge(n,2);
    n2        = edge(n,3);
    xa        = x(n1,0);
    ya        = x(n1,1);
    xb        = x(n2,0);
    yb        = x(n2,1);
    area(n,0) = ya-yb;
    area(n,1) = xb-xa;
  }
  for (int n=nEdge-nEdgeBd; n<nEdge; n++){ //impinging face areas
    p1        = edge(n,0);
    p2        = edge(n,1);
    n1        = edge(n,2);
    xa        = x(n1,0);
    ya        = x(n1,1);
    xb        = .5*(x(p1,0)+x(p2,0));
    yb        = .5*(x(p1,1)+x(p2,1));
    area(n,0) = ya-yb;
    area(n,1) = xb-xa;
  }
  for (int n=0; n<nEdge; n++)
    for (int j=0; j<2; j++) area(n,j) /= 3.;
  int m=0;
  for (int n=nEdge-nEdgeBd; n<nEdge; n++){ //boundary face areas
    p1            = edge(n,0);
    p2            = edge(n,1);
    areaBd(m  ,0) = .5*(x(p2,1)-x(p1,1));
    areaBd(m++,1) = .5*(x(p1,0)-x(p2,0));
  }

  for (int n=0; n<nEdge; n++)
    for (int j=0; j<2; j++) area(n,j) *= .5;
  for (int n=0; n<nEdgeBd; n++)
    for (int j=0; j<2; j++) areaBd(n,j) *= .5;


  Array2D<double> sumA(nNode,2); //test for zero sum areas
  sumA.set(0.);
  for (int n=0; n<nEdge; n++){
    p1         = edge(n,0);
    p2         = edge(n,1);
    sumA(p1,0) = sumA(p1,0)+area(n,0);
    sumA(p1,1) = sumA(p1,1)+area(n,1);
    sumA(p2,0) = sumA(p2,0)-area(n,0);
    sumA(p2,1) = sumA(p2,1)-area(n,1);
  }
  m = 0;
  for (int n=nEdge-nEdgeBd; n<nEdge; n++){
    p1         = edge(n,0);
    p2         = edge(n,1);
    sumA(p1,0) = sumA(p1,0)+areaBd(m  ,0);
    sumA(p1,1) = sumA(p1,1)+areaBd(m  ,1);
    sumA(p2,0) = sumA(p2,0)+areaBd(m  ,0);
    sumA(p2,1) = sumA(p2,1)+areaBd(m++,1);
  }
  xa = 0.;
  ya = 0.;
  for (int n=0; n<nNode; n++){
    xa = max(xa,fabs(sumA(n,0)));
    ya = max(ya,fabs(sumA(n,1)));
  }
  cout << "\nmax face area sum (" << xa << ", " << ya << ")" << endl;
  sumA.deallocate();
}
