#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::volume()
{
  v.set(0.);
  int p1,p2,p3;
  double xa,ya,xb,yb,xc,yc,vv;
  for (int n=0; n<nTri; n++){
    p1     = tri(n,0);
    p2     = tri(n,1);
    p3     = tri(n,2);
    xa     = x(p1,0);
    ya     = x(p1,1);
    xb     = x(p2,0);
    yb     = x(p2,1);
    xc     = x(p3,0);
    yc     = x(p3,1);
    vv     = .5*((xb-xa)*(yc-ya)-(yb-ya)*(xc-xa));
    v(p1) += vv;
    v(p2) += vv;
    v(p3) += vv;
  }
  for (int n=0; n<nNode; n++) v(n) /= 6.;

  vv = 0.;
  for (int n=0; n<nNode; n++) vv += v(n);
  cout << "\nFine level total volume: " << vv << endl;
}
