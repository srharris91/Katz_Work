#include "StrandBlockSolver.h"
#include "edgeExtract.h"
#include "formCsp.h"
#include "faceArea.h"
#include "cellCoord.h"
#include "volume.h"


void StrandBlockSolver::setStrandBlockData(StrandBlock& strandBlock)
{
  //obtain data from this strandBlock
  pid                  = strandBlock.getPid();
  nFaces               = strandBlock.getNFaces();
  nGfaces              = strandBlock.getNGfaces();
  nNodes               = strandBlock.getNNodes();
  nGnodes              = strandBlock.getNGnodes();
  nBedges              = strandBlock.getNBedges();
  nPedges              = strandBlock.getNPedges();
  nPstr                = strandBlock.getNPstr();
  nFringe              = strandBlock.getNFringe();
  nSharp               = strandBlock.getNSharp();
  const int*    faceT  = strandBlock.getFace();
  const int*    fTagT  = strandBlock.getFTag();
  const int*    sFlagT = strandBlock.getSFlag();
  const int*    bEdgeT = strandBlock.getBEdge();
  const int*    pEdgeT = strandBlock.getPEdge();
  const int*    bTagT  = strandBlock.getBTag();
  const int*    fClipT = strandBlock.getFClip();
  const double* xSrfT  = strandBlock.getXSrf();
  const double* xStrT  = strandBlock.getXStr();
  const double* pVT    = strandBlock.getPV();


  // allocate local solver data pertaining to connectivity
  // and geometry
  int i,j,k,l;
  i = nFaces+nBedges;
  l = nFaces+nBedges+1;
  j = nPstr+1;
  k = nPstr+2;
  nEdges =(2*(nFaces-nGfaces)+nBedges+nPedges)/2;
  face.allocate(2,nFaces);
  fTag.allocate(nFaces);
  bTag.allocate(nBedges);
  fClip.allocate(i);
  sFlag.allocate(nNodes);
  edge.allocate(2,nEdges);
  edgp.allocate(2,nEdges);
  edgn.allocate(nEdges);
  ncsc.allocate(l);
  ncsp.allocate(nNodes);
  gsMap.allocate(i);
  x.allocate(ndim,j,nNodes);
  xc.allocate(ndim,k,i);
  v.allocate(k,i);
  facu.allocate(ndim,j,nFaces);
  facs.allocate(ndim,j,nEdges);
  xStr.allocate(j);
  csp = new int*[nNodes];
  lsp = new double*[2*j*nNodes];


  // fill in what I can that is easy
  for (int n=0; n<nFaces ; n++){
    face(0,n) = faceT[2*n  ]-1; //0-based
    face(1,n) = faceT[2*n+1]-1; //0-based
  }
  for (int n=0; n<  nFaces ; n++) fTag(n) = fTagT[n]-1;//0-based
  for (int n=0; n<  nBedges; n++) bTag(n) = bTagT[n]-1;//0-based
  for (int n=0; n<  nFaces ; n++) fClip(n) = fClipT[n];
  for (int n=0; n<nFaces+nBedges; n++) fClip(n) = nPstr-nFringe;
  if (standAlone == 1) for (int n=0; n<nFaces+nBedges; n++)
			 fClip(n) = nPstr;
  for (int n=0; n<  nNodes ; n++) sFlag(n) = sFlagT[n];
  for (int n=0; n<  nPstr+1; n++) xStr(n) = xStrT[n];
  double xn,yn,nx,ny;
  for (int n=0; n<nNodes; n++){
    xn  = xSrfT[2*n  ];
    yn  = xSrfT[2*n+1];
    nx  = pVT  [2*n  ];
    ny  = pVT  [2*n+1];
    for (int j=0; j<nPstr+1; j++){
      x(0,j,n) = xn+xStr(j)*nx;
      x(1,j,n) = yn+xStr(j)*ny;
    }
  }

  cout << "\nFine level clipping:" << endl;
  for (int n=0; n<nFaces+nBedges; n++)
    cout << n << "\t" << fClip(n) << endl;


  // extract edges
  edgeextract_(pid,
	       nFaces,
	       nGfaces,
	       nNodes,
	       nGnodes,
	       nEdges,
	       nBedges,
	       nPedges,
	       &face(0,0),
	       bEdgeT,
	       pEdgeT,
	       &bTag(0),
	       &edge(0,0),
	       &edgp(0,0),
	       &edgn(0));


  // perturb nodes for verification purposes only
  if (perturb != 0) perturbNodes();


  // form cells surrounding points
  int  ncsp1 = 3*nNodes; // hopefully enough
  int* csp1  = new int[ncsp1   ];
  int* csp2  = new int[nNodes+1];
  formcsp_(pid,
	   nFaces,
	   nGfaces,
	   nNodes,
	   nGnodes,
	   ncsp1,
	   &face(0,0),
	   csp1,
	   csp2);
  for (int n=0; n<nNodes; n++){
    ncsp(n) = csp2[n+1]-csp2[n];
    csp[n]  = new int[ncsp(n)];
    for (int m=0; m<ncsp(n); m++) csp[n][m] = csp1[csp2[n]+m];
  }
  delete [] csp1;
  delete [] csp2;


  // set cells surrounding cells for implicit neighbors
  formCsc();


  // set cell ordering for GS procedure
  order();


  // face areas
  facearea_(pid,
	    ndim,
	    nFaces,
	    nGfaces,
	    nNodes,
	    nGnodes,
	    nEdges,
	    nBedges,
	    nPstr,
	    &face(0,0),
	    &edge(0,0),
	    &edgn(0),
	    &x(0,0,0),
	    &facs(0,0,0),
	    &facu(0,0,0));


  // cell-center coordinates
  cellcoord_(pid,
	     ndim,
	     nFaces,
	     nGfaces,
	     nNodes,
	     nGnodes,
	     nEdges,
	     nBedges,
	     nPstr,
	     &face(0,0),
	     &edge(0,0),
	     &edgn(0),
	     &x(0,0,0),
	     &facs(0,0,0),
	     &facu(0,0,0),
	     &xc(0,0,0));


  // cell volumes
  volume_(pid,
	  ndim,
	  nFaces,
	  nGfaces,
	  nNodes,
	  nGnodes,
	  nEdges,
	  nBedges,
	  nPstr,
	  &face(0,0),
	  &edge(0,0),
	  &edgn(0),
	  &fClip(0),
	  &x(0,0,0),
	  &facs(0,0,0),
	  &facu(0,0,0),
	  &xc(0,0,0),
	  &v(0,0));


  // least squares interpolation coefficients
  int ii,mm;
  for (int n=0; n<nNodes; n++){
    mm = ncsp(n);
    for (int j=0; j<nPstr+1; j++){
      for (int k=0; k<2; k++){
	indlsp(k,j,n,ii);
	lsp[ii] = new double[mm];
      }
    }
  }
  if      (nodeVal == 1) lspVol();
  else if (nodeVal == 2) lspLS();
  else if (nodeVal == 3) lspMap();
  else{
    cout << "\nvalue of nodeVal not recognized" << endl;
    exit(0);
  }

  // test the ability of the interplation to compute constant and linear
  // functions
  int jj,nn;
  double sum,sumx,sumy,sumM,sumxM,sumyM;
  sumM  = 0.;
  sumxM = 0.;
  sumyM = 0.;
  for (int n=0; n<nNodes-nGnodes; n++){
    mm = ncsp(n);
  for (int j=0; j<nPstr+1; j++){
    if      (j == 0    ) jj = 1;
    else if (j == nPstr) jj = nPstr-1;
    else                 jj = j;
    sum  = 0.;
    sumx = 0.;
    sumy = 0.;
  for (int k=0; k<2; k++){
    indlsp(k,j,n,ii);
  for (int m=0; m<mm; m++){
    nn    = csp[n][m];
    sum  += lsp[ii][m];
    sumx += lsp[ii][m]*xc(0,jj,nn);
    sumy += lsp[ii][m]*xc(1,jj,nn);
  }
  jj++;
  }
    sum   = fabs(sum -1.      );
    sumx  = fabs(sumx-x(0,j,n));
    sumy  = fabs(sumy-x(1,j,n));
    if (sum  > sumM ) sumM  = sum;
    if (sumx > sumxM) sumxM = sumx;
    if (sumy > sumyM) sumyM = sumy;
  }
  }
  cout << "\n";
  cout << "Constant function lsp max: " << sumM  << endl;
  cout << "Linear-x function lsp max: " << sumxM << endl;
  cout << "Linear-y function lsp max: " << sumyM << endl;
  cout << "\n";

  
  // check that the number of boundary patches equals the input file
  i = 0;
  for (int n=0; n<nFaces ; n++) if (fTag(n) > i) i = fTag(n);
  for (int n=0; n<nBedges; n++) if (bTag(n) > i) i = bTag(n);
  i++; //strand ends
  i++; //inputs are 1-based
  if (i != nBpatches){
    cout << "\nNumber of boundary patches specified in input file does ";
    cout << "not correspond to the grid file. Terminating." <<" nBpatches = "<<nBpatches<<" i = "<<i<< endl;
    exit(0);
  }
}
