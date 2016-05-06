#include "StrandBlockSolver.h"
#include "coarsenSBS.h"


void StrandBlockSolver::coarsen(const int& mglevel,
				StrandBlockSolver& parentSbs)
{
  // set data from the parent block solver
  pid      = parentSbs.getPid();
  nq       = parentSbs.getNq();
  nqa      = parentSbs.getNqa();
  ndim     = parentSbs.getNdim();
  nBpatches = parentSbs.getNBpatches();
  inviscid = parentSbs.getInviscid();
  viscous  = parentSbs.getViscous();
  source   = parentSbs.getSource();
  dissipation = parentSbs.getDissipation();
  sys      = parentSbs.getSystem();


  // set pointers to parent block solver data
  limFlagP = parentSbs.getLimFlag();
  nodalQFlagP = parentSbs.getNodalQFlag();
  nodalQaFlagP = parentSbs.getNodalQaFlag();
  gradQFlagP = parentSbs.getGradQFlag();
  gradQaFlagP = parentSbs.getGradQaFlag();
  nFacesP  = parentSbs.getNFaces();
  nGfacesP = parentSbs.getNGfaces();
  nEdgesP  = parentSbs.getNEdges();
  nBedgesP = parentSbs.getNBedges();
  nPedgesP = parentSbs.getNPedges();
  nPstrP   = parentSbs.getNPstr();
  nFringeP = parentSbs.getNFringe();
  edgeP    = parentSbs.getEdgeArray();
  bTagP    = parentSbs.getBTagArray();
  fTagP    = parentSbs.getFTagArray();
  fClipP   = parentSbs.getFClipArray();
  facuP    = parentSbs.getFacuArray();
  facsP    = parentSbs.getFacsArray();
  xvuP     = parentSbs.getXvuArray();
  xvsP     = parentSbs.getXvsArray();
  qP       = parentSbs.getQArray();
  qaP      = parentSbs.getQaArray();
  rP       = parentSbs.getRArray();
  vP       = parentSbs.getVArray();


  // allocate data necessary for the coarsening procedure (this will involve
  // overdimensioning some arrays at first)
  int i,j,k,l;
  j = nFacesP+nBedgesP;
  k = nPstrP+1;
  l = nPstrP+2;
  f2cc.allocate(j);
  f2ce.allocate(nEdgesP);
  f2cs.allocate(l);
  Array1D<int> fTagT(nFacesP);
  Array1D<int> bTagT(nBedgesP);
  Array1D<int> fClipT(j);
  Array2D<int> edgeT(ndim,nEdgesP);
  int nFacesT;
  int nGfacesT;
  int nPstrT;
  int nBedgesT;
  int nPedgesT;
  int nEdgesT;

  coarsensbs_(standAlone,
	      mglevel,
	      nFacesT,
	      nGfacesT,
	      nPstrT,
	      nBedgesT,
	      nPedgesT,
	      nEdgesT,
	      nFacesP,
	      nGfacesP,
	      nEdgesP,
	      nBedgesP,
	      nPedgesP,
	      nPstrP,
	      nFringeP,
	      ndim,
	      &(*edgeP)(0,0),
	      &(*bTagP)(0),
	      &(*fTagP)(0),
	      &(*fClipP)(0),
	      &f2cc(0),
	      &f2ce(0),
	      &f2cs(0),
	      &fTagT(0),
	      &bTagT(0),
	      &fClipT(0),
	      &edgeT(0,0));

  nFaces  = nFacesT;
  nGfaces = nGfacesT;
  nPstr   = nPstrT;
  nBedges = nBedgesT;
  nPedges = nPedgesT;
  nEdges  = nEdgesT;
  nFringe = 1;

  j = nFaces+nBedges;
  i = nFaces+nBedges+1;
  k = nPstr+1;
  l = nPstr+2;
  fTag.allocate(nFaces);
  bTag.allocate(nBedges);
  fClip.allocate(j);
  edge.allocate(ndim,nEdges);

  for (int n=0; n<nFaces; n++) fTag(n) = fTagT(n);
  for (int n=0; n<nBedges; n++) bTag(n) = bTagT(n);
  for (int n=0; n<nFaces+nBedges; n++) fClip(n) = fClipT(n);
  for (int n=0; n<nEdges; n++){
    edge(0,n) = edgeT(0,n);
    edge(1,n) = edgeT(1,n);
  }

  fTagT.deallocate();
  bTagT.deallocate();
  fClipT.deallocate();
  edgeT.deallocate();

  gsMap.allocate(j);
  v.allocate(l,j);
  facu.allocate(ndim,k,nFaces);
  facs.allocate(ndim,k,nEdges);
  ncsc.allocate(i);


  // set cells surrounding cells for implicit neighbors
  formCsc();


  // set cell ordering for GS procedure
  order();


  // compute coarse level face areas and volumes
  coarseMetrics();


  // check volume sum and face area sum
  cout << "\nChecking grid level: " << mglevel << endl;
  int jp;
  double a=0.;
  for (int n=0; n<nFaces-nGfaces; n++)
    for (int j=1; j<nPstr+1; j++) a += v(j,n);
  cout << "Coarse level total volume: " << a << endl;

  j = nFaces+nBedges;
  k = nPstr+1;
  l = nPstr+2;
  Array3D<double> A(ndim,l,j);
  int c1,c2;
  double xm=0.,ym=0.;
  for (int n=0; n<nFaces+nBedges; n++)
    for (int j=0; j<nPstr+2; j++){
      A(0,j,n) = 0.;
      A(1,j,n) = 0.;
    }
  for (int n=0; n<nEdges; n++){
    c1 = edge(0,n);
    c2 = edge(1,n);
  for (int j=1; j<nPstr+1; j++){
    A(0,j,c1) += facs(0,j,n);
    A(1,j,c1) += facs(1,j,n);
    A(0,j,c2) -= facs(0,j,n);
    A(1,j,c2) -= facs(1,j,n);
  }}
  for (int n=0; n<nFaces-nGfaces; n++){
  for (int j=0; j<nPstr+1; j++){
    jp         = j+1;
    A(0,j ,n) += facu(0,j,n);
    A(1,j ,n) += facu(1,j,n);
    A(0,jp,n) -= facu(0,j,n);
    A(1,jp,n) -= facu(1,j,n);
  }}
  for (int n=0; n<nFaces-nGfaces; n++)
    for (int j=1; j<fClip(n)+1; j++){
      if (fabs(A(0,j,n)) > xm) xm = fabs(A(0,j,n));
      if (fabs(A(1,j,n)) > ym) ym = fabs(A(1,j,n));
    }
  cout << "Max x-face area sum: " << xm << endl;
  cout << "Max y-face area sum: " << ym << endl;
}
