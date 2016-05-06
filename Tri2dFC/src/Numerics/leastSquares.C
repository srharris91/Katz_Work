#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::leastSquares()
{
  //form psp arrays
  int i,j;
  psp2.allocate(nNode+1);
  psp2.set(0);
  for (int n=0; n<nNode; n++) psp2(n+1) = psp2(n+1)+1;
  for (int n=0; n<nEdge; n++){
    j = edge(n,0)+1;
    psp2(j)++;
    j = edge(n,1)+1;
    psp2(j)++;
  }
  for (int n=1; n<nNode+1; n++) psp2(n) += psp2(n-1);

  npsp1 = psp2(nNode);
  psp1.allocate(npsp1);
  for (int n=0; n<nNode; n++){
    i       = psp2(n);
    psp1(i) = n;
    psp2(n)++;
  }
  for (int n=0; n<nEdge; n++){
    j       = edge(n,0);
    i       = psp2(j);
    psp1(i) = edge(n,1);
    psp2(j)++;
    j       = edge(n,1);
    i       = psp2(j);
    psp1(i) = edge(n,0);
    psp2(j)++;
  }
  for (int n=nNode; n>0; n--) psp2(n) = psp2(n-1);
  psp2(0) = 0;

  /*
  for (int n=0; n<nNode; n++){
    cout << "\nNode: " << n << " " << psp2(n+1)-psp2(n) << endl;
    for (int i=psp2(n); i<psp2(n+1); i++) cout << psp1(i) << endl;
  }
  exit(0);
  */


  //form linear lsq coefficients and store in gx
  double xi,yi,xn,yn,dsmax,w,sxx,sxy,syy,dx,dy;
  gx.allocate(npsp1,2);
  gx.set(0.);
  for(int n=0; n<nNode; n++){
    xn = x(n,0);
    yn = x(n,1);

    dsmax = 0.; // Find largest distance in stencil
    for (int m=psp2(n)+1; m<psp2(n+1); m++){
      i     = psp1(m);
      xi    = x(i,0);
      yi    = x(i,1);
      dx    = xi-xn;
      dy    = yi-yn;
      dsmax = max(dsmax,dx*dx+dy*dy);
    }
    dsmax = 1./sqrt(dsmax);

    sxx = 0.; // Form least squares matrix
    sxy = 0.;
    syy = 0.;
    for (int m=psp2(n)+1; m<psp2(n+1); m++){
      i    = psp1(m);
      xi   = x(i,0);
      yi   = x(i,1);
      dx   =(xi-xn);//*dsmax;
      dy   =(yi-yn);//*dsmax;
      w    = 1.;///(dx*dx+dy*dy);
      sxx += w*dx*dx;
      sxy += w*dx*dy;
      syy += w*dy*dy;
    }

    // Compute and store least squares coefficients
    w    = 1./(sxx*syy-sxy*sxy);
    sxx *= w;
    sxy *= w;
    syy *= w;
    for (int m=psp2(n)+1; m<psp2(n+1); m++){
      i       = psp1(m);
      xi      = x(i,0);
      yi      = x(i,1);
      dx      =(xi-xn);//*dsmax;
      dy      =(yi-yn);//*dsmax;
      w       = 1.;///(dx*dx+dy*dy);
      gx(m,0) =( syy*w*dx-sxy*w*dy);//*dsmax;
      gx(m,1) =(-sxy*w*dx+sxx*w*dy);//*dsmax;
    }
    dx = 0.;
    dy = 0.;
    for (int m=psp2(n)+1; m<psp2(n+1); m++){
      dx += gx(m,0);
      dy += gx(m,1);
    }
    gx(psp2(n),0) =-dx;
    gx(psp2(n),1) =-dy;
  }


  // check coefficients for errors
  Array1D<double> sumax,sumay,sumxx,sumxy,sumyx,sumyy;
  sumax.allocate(nNode);
  sumay.allocate(nNode);
  sumxx.allocate(nNode);
  sumxy.allocate(nNode);
  sumyx.allocate(nNode);
  sumyy.allocate(nNode);
  sumax.set(0.);
  sumay.set(0.);
  sumxx.set(0.);
  sumxy.set(0.);
  sumyx.set(0.);
  sumyy.set(0.);

  for(int n=0; n<nNode; n++)
    for (int m=psp2(n); m<psp2(n+1); m++){
      i         = psp1(m);
      sumax(n) += gx(m,0);
      sumay(n) += gx(m,1);
      sumxx(n) += gx(m,0)*x(i,0);
      sumxy(n) += gx(m,0)*x(i,1);
      sumyx(n) += gx(m,1)*x(i,0);
      sumyy(n) += gx(m,1)*x(i,1);
    }

  cout << "\n" << endl;

  dx = 10.;
  dy =-10.;
  for(int n=0; n<nNode; n++){
    dx = min(dx,sumax(n));
    dy = max(dy,sumax(n));
  }
  cout << "max of ls sumax: " << dy << endl;
  cout << "min of ls sumax: " << dx << endl;

  dx = 10.;
  dy =-10.;
  for(int n=0; n<nNode; n++){
    dx = min(dx,sumay(n));
    dy = max(dy,sumay(n));
  }
  cout << "max of ls sumay: " << dy << endl;
  cout << "min of ls sumay: " << dx << endl;

  dx = 10.;
  dy =-10.;
  for(int n=0; n<nNode; n++){
    dx = min(dx,sumxx(n));
    dy = max(dy,sumxx(n));
  }
  cout << "max of ls sumxx: " << dy << endl;
  cout << "min of ls sumxx: " << dx << endl;

  dx = 10.;
  dy =-10.;
  for(int n=0; n<nNode; n++){
    dx = min(dx,sumxy(n));
    dy = max(dy,sumxy(n));
  }
  cout << "max of ls sumxy: " << dy << endl;
  cout << "min of ls sumxy: " << dx << endl;

  dx = 10.;
  dy =-10.;
  for(int n=0; n<nNode; n++){
    dx = min(dx,sumyx(n));
    dy = max(dy,sumyx(n));
  }
  cout << "max of ls sumyx: " << dy << endl;
  cout << "min of ls sumyx: " << dx << endl;

  dx = 10.;
  dy =-10.;
  for(int n=0; n<nNode; n++){
    dx = min(dx,sumyy(n));
    dy = max(dy,sumyy(n));
  }
  cout << "max of ls sumyy: " << dy << endl;
  cout << "min of ls sumyy: " << dx << endl;

  cout << "\n" << endl;

  sumax.deallocate();
  sumay.deallocate();
  sumxx.deallocate();
  sumxy.deallocate();
  sumyx.deallocate();
  sumyy.deallocate();
}
