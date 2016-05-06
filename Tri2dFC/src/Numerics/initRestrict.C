#include "Tri2dFCBlockSolver.h"
#include "lagrangePoly.h"


void Tri2dFCBlockSolver::initRestrict(int levelF0,
				      int nNodeF0,
				      int nNodeBdF0,
				      int* elemF0,
				      int* nodeBdF,
				      double* xF,
				      double* qF0,
				      double* rF0)
{
  // set the pointers to the fine data
  levelF   = levelF0;
  nNodeF   = nNodeF0;
  nNodeBdF = nNodeBdF0;
  elemF    = elemF0;
  qF       = qF0;
  rF       = rF0;


  // for each coarse node, find the fine level element which contains it,
  // as well as its (r,s) index in the standard element
  nfe.allocate(nNode,2);
  nfe.set(-1);
  for (int n=0; n<nElem; n++)
    for (int j=0; j<nne; j++)
      if (nfe(elem(n,j),0) == -1){
	nfe(elem(n,j),0) = n;
	nfe(elem(n,j),1) = j;
      }


  // set up solution restriction coefficients via interpolation
  // lqFC(i,j) = l_j(r_i) (a row is all fine level Lagrange polynomials
  // evaluated at a single coarse level solution point i)
  int orderE=3-level; //order of local elements on this level
  int orderEF=3-levelF; //order of local elements on fine level
  bool test=false;
  nneF =(orderEF+2)*(orderEF+1)/2; //number of nodes per element on fine level
  Array2D<double> rs(nne,3),rsF(nneF,3),lcF(nneF,nneF);
  solutionPoints(orderE,
		 spacing,
		 &rs(0,0));
  solutionPoints(orderEF,
		 spacing,
		 &rsF(0,0));
  lagrangePoly(test,
	       orderEF,
	       &rsF(0,0),
	       &lcF(0,0));
  lqFC.allocate(nne,nneF);
  lqFC.set(0.);
  int j;
  double ri,si;
  for (int n=0; n<nneF; n++) // nth fine level Lagrange polynomial
    for (int i=0; i<nne; i++){ // ith coarse level solution point
      j = 0;
      for (int k=0; k<=orderEF; k++)
	for (int l=0; l<=orderEF-k; l++){
	  ri = rs(i,0);
	  si = rs(i,1);
	  lqFC(i,n) += pow(ri,k)*pow(si,l)*lcF(n,j++);
	}
    }

  rs.deallocate();
  rsF.deallocate();
  lcF.deallocate();

  /*
  for (int i=0; i<nne; i++){
    for (int j=0; j<nneF; j++) cout << lqFC(i,j) << " ";
    cout << endl;
  }
  exit(0);
  */


  // set up residual restriction coefficients
  // first flag the boundary nodes with their boundary tags
  Array1D<int> flag(nNode);
  flag.set(-1);
  int m=0;
  for (int n=nNode-nNodeBd; n<nNode; n++) flag(n) = nodeBd(m++);
  Array1D<int> flagF(nNodeF);
  flagF.set(-1);
  m = 0;
  for (int n=nNodeF-nNodeBdF; n<nNodeF; n++) flagF(n) = nodeBdF[m++];

  // then allocate and initialize the transfer arrays
  nfn  = new int    [nNodeF];
  nfn1 = new int*   [nNodeF];
  nfn2 = new double*[nNodeF];
  for (int n=0; n<nNodeF; n++) nfn [n] = 0;
  for (int n=0; n<nNodeF; n++) nfn1[n] = NULL;
  for (int n=0; n<nNodeF; n++) nfn2[n] = NULL;
  int npF,n0C,n1C,n2C,triL[3];
  double xpF,ypF,x0C,y0C,x1C,y1C,x2C,y2C,A0,A1,A2,A,eps=1.e-14;
  for (int n=0; n<nElem; n++){
    for (int jF=0; jF<nneF; jF++){
      npF = elemF[n*nneF+jF];
      if (nfn[npF] == 0){
	if (level == 1){
	  if      (jF == 0){
	    triL[0] = 0;
	    triL[1] = 3;
	    triL[2] = 5;
	  }
	  else if (jF == 1){
	    triL[0] = 3;
	    triL[1] = 1;
	    triL[2] = 4;
	  }
	  else if (jF == 2){
	    triL[0] = 5;
	    triL[1] = 4;
	    triL[2] = 2;
	  }
	  else if (jF == 3){
	    triL[0] = 0;
	    triL[1] = 3;
	    triL[2] = 5;
	  }
	  else if (jF == 4){
	    triL[0] = 3;
	    triL[1] = 1;
	    triL[2] = 4;
	  }
	  else if (jF == 5){
	    triL[0] = 3;
	    triL[1] = 1;
	    triL[2] = 4;
	  }
	  else if (jF == 6){
	    triL[0] = 5;
	    triL[1] = 4;
	    triL[2] = 2;
	  }
	  else if (jF == 7){
	    triL[0] = 5;
	    triL[1] = 4;
	    triL[2] = 2;
	  }
	  else if (jF == 8){
	    triL[0] = 0;
	    triL[1] = 3;
	    triL[2] = 5;
	  }
	  else if (jF == 9){
	    triL[0] = 3;
	    triL[1] = 4;
	    triL[2] = 5;
	  }}
	else if (level == 2){
	  triL[0] = 0;
	  triL[1] = 1;
	  triL[2] = 2;
	}
	n0C = elem(n,triL[0]);
	n1C = elem(n,triL[1]);
	n2C = elem(n,triL[2]);
	xpF = xF[npF*2  ];
	ypF = xF[npF*2+1];
	x0C = x(n0C,0)-xpF;
	y0C = x(n0C,1)-ypF;
	x1C = x(n1C,0)-xpF;
	y1C = x(n1C,1)-ypF;
	x2C = x(n2C,0)-xpF;
	y2C = x(n2C,1)-ypF;
	A0  = .5*fabs(x1C*y2C-x2C*y1C);
	A1  = .5*fabs(x2C*y0C-x0C*y2C);
	A2  = .5*fabs(x0C*y1C-x1C*y0C);
	A   = A0+A1+A2;
	A0 /= A;
	A1 /= A;
	A2 /= A;
	if (A0 > eps && flag(n0C) == flagF(npF)) nfn[npF]++;
	if (A1 > eps && flag(n1C) == flagF(npF)) nfn[npF]++;
	if (A2 > eps && flag(n2C) == flagF(npF)) nfn[npF]++;
	if (nfn[npF] > 0){
	  nfn1[npF] = new int   [nfn[npF]];
	  nfn2[npF] = new double[nfn[npF]];
	  nfn [npF] = 0;
	  if (A0 > eps && flag(n0C) == flagF(npF)){
	    nfn1[npF][nfn[npF]  ] = n0C;
	    nfn2[npF][nfn[npF]++] = A0;
	  }
	  if (A1 > eps && flag(n1C) == flagF(npF)){
	    nfn1[npF][nfn[npF]  ] = n1C;
	    nfn2[npF][nfn[npF]++] = A1;
	  }
	  if (A2 > eps && flag(n2C) == flagF(npF)){
	    nfn1[npF][nfn[npF]  ] = n2C;
	    nfn2[npF][nfn[npF]++] = A2;
	  }}}}}

  // normalize boundary transfer coefficients
  Array1D<double> sumC(nNode);
  sumC.set(0.);
  int i;
  double a;
  for (int n=nNodeF-nNodeBdF; n<nNodeF; n++)
    for (int j=0; j<nfn[n]; j++){
      i        = nfn1[n][j];
      a        = nfn2[n][j];
      sumC(i) += a;
    }
  for (int n=nNodeF-nNodeBdF; n<nNodeF; n++)
    for (int j=0; j<nfn[n]; j++){
      i           = nfn1[n][j];
      nfn2[n][j] /= sumC(i);
    }

  flag.deallocate();
  flagF.deallocate();
  sumC.deallocate();

  /*
  for (int n=0; n<nNodeF; n++){
    cout << n << " ";
    for (int j=0; j<nfn[n]; j++){
      i       = nfn1[n][j];
      a       = nfn2[n][j];
      cout << i << " " << a << " ";
    }
    cout << endl;
  }
  */
}
