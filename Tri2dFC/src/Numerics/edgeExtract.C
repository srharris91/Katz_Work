#include "Tri2dFCBlockSolver.h"


void Tri2dFCBlockSolver::edgeExtract()
{
  // form cells surrounding points
  int* csp2 = new int[nNode];
  for (int n=0; n<nNode; n++) csp2[n] = 0;
  for (int n=0; n<nTri; n++)
    for (int k=0; k<3; k++) csp2[tri(n,k)]++;

  int** csp1 = new int*[nNode];
  for (int n=0; n<nNode; n++){
    if (csp2[n] > 0) csp1[n] = new int[csp2[n]];
    else csp1[n] = NULL;
  }

  int j;
  for (int n=0; n<nNode; n++) csp2[n] = 0;
  for (int n=0; n<nTri; n++)
    for (int k=0; k<3; k++){
      j = csp2[tri(n,k)];
      csp1[tri(n,k)][j] = n;
      csp2[tri(n,k)]++;
    }

  // for (int n=0; n<nNode; n++){
  //   cout << n << " ";
  //   for (int k=0; k<csp2[n]; k++)
  //     cout << csp1[n][k] << " ";
  //   cout << endl;
  // }
  // exit(0);


  // extract edges
  // edges have the following ordering:
  // edge(n,0) - first node on edge
  // edge(n,1) - second node on edge
  // edge(n,2) - CCW opposing node on edge
  // edge(n,3) - CW opposing node on edge (or boundary tag number)
  // edge(n,4) - edge extension node for node 1
  // edge(n,5) - edge extension node for node 2
  nEdge =(3*nTri+nEdgeBd)/2;
  edge.allocate(nEdge,6);
  edge.set(-1);
  Array2D<int> nflag(nNode,2);
  nflag.set(-1);
  int i,l,k1,mm,n1,n2,ne=0;
  for (int n=0; n<nNode; n++){ //for each node
    for (int m=0; m<csp2[n]; m++){ //look at the cells surrounding the node
      j = csp1[n][m];
      for (int k=0; k<3; k++) //find the position of n in tri
	if (n == tri(j,k)){
	  k1 = k;
	  break;
	}
      for (int k=0; k<3; k++){
	l = tri(j,k);
	if (n < l){
	  if ((k1 == 0 && k == 1) ||
	      (k1 == 1 && k == 2) ||
	      (k1 == 2 && k == 0)) n1 = 0; //n to l is CCW
	  else n1 = 1; //n to l is CW
	  for (int kk=0; kk<3; kk++){
            mm = tri(j,kk); //opposing node
	    if (mm != n && mm != l) break;
	  }
	  if (nflag(l,0) != n){ //if haven't added this edge yet
	    if (n1 == 0){//n to l is CCW
	      edge(ne,0) = n;
	      edge(ne,1) = l;
	      edge(ne,2) = mm;
	      edge(ne,4) = j;
	    }
	    else{//n to l is CW
	      edge(ne,1) = n;
	      edge(ne,0) = l;
	      edge(ne,2) = mm;
	      edge(ne,4) = j;
	    }
	    nflag(l,0) = n;
	    nflag(l,1) = ne;
	    ne++;
	  }
	  else{//have added this edge, just add the opposing node
	    i         = nflag(l,1);
	    edge(i,3) = mm;
	    edge(i,5) = j;
	  }}}}}

  if (ne != nEdge){
    cout << "\n*** number of edges counted incorrectly in connectivity.C ***"
	 << endl;
    exit(0);
  }

  // for (int n=0; n<nEdge; n++)
  //   cout << n << " "
  // 	 << edge(n,0) << " " << edge(n,1) << " "
  // 	 << edge(n,2) << " " << edge(n,3) << " "
  // 	 << edge(n,4) << " " << edge(n,5) << endl;
  // exit(0);


  // put boundary edges last
  Array2D<int> edgeT(nEdge,6);
  ne = 0;
  for (int n=0; n<nEdge; n++)
    if (edge(n,3) >= 0){
      edgeT(ne  ,0) = edge(n,0);
      edgeT(ne  ,1) = edge(n,1);
      edgeT(ne  ,2) = edge(n,2);
      edgeT(ne  ,3) = edge(n,3);
      edgeT(ne  ,4) = edge(n,4);
      edgeT(ne++,5) = edge(n,5);
    }

  if (ne != nEdge-nEdgeBd){
    cout << "\n*** number of b. edges counted incorrectly in connectivity.C ***"
	 << endl;
    exit(0);
  }

  for (int n=0; n<nEdge; n++)
    if (edge(n,3) < 0){
      edgeT(ne  ,0) = edge(n,0);
      edgeT(ne  ,1) = edge(n,1);
      edgeT(ne  ,2) = edge(n,2);
      edgeT(ne  ,3) = edge(n,3);
      edgeT(ne  ,4) = edge(n,4);
      edgeT(ne++,5) = edge(n,5);
    }

  if (ne != nEdge){
    cout << "\n*** number of edges counted incorrectly in connectivity.C ***"
	 << endl;
    exit(0);
  }

  for (int n=0; n<nEdge; n++){
    edge(n,0) = edgeT(n,0);
    edge(n,1) = edgeT(n,1);
    edge(n,2) = edgeT(n,2);
    edge(n,3) = edgeT(n,3);
    edge(n,4) = edgeT(n,4);
    edge(n,5) = edgeT(n,5);
  }


  // flag boundary nodes with the two boundary edge numbers that share it
  nflag.set(-1);
  for (int n=nEdge-nEdgeBd; n<nEdge; n++){
    n1 = edge(n,0);
    n2 = edge(n,1);
    if (nflag(n1,0) == -1) nflag(n1,0) = n;
    else nflag(n1,1) = n;
    if (nflag(n2,0) == -1) nflag(n2,0) = n;
    else nflag(n2,1) = n;
  }


  // for each boundary edge fill in the boundary tag
  for (int n=0; n<nEdgeBd; n++){
    n1 = edgeBd(n,0);
    n2 = edgeBd(n,1);
    if      (nflag(n1,0) == nflag(n2,0) ||
	     nflag(n1,0) == nflag(n2,1)) k1 = nflag(n1,0);
    else k1 = nflag(n1,1);
    if (edge(k1,3) != -1){
      cout << "\n*** error identifying boundary tags in connectivity.C ***"
	   << endl;
      exit(0);
    }
    else edge(k1,3) = edgeBd(n,2);
  }


  // find extension triangles
  int k2;
  double dx,dy,ds,dx1,dy1,dx2,dy2,cp1,cp2;
  for (int n=nEdge-nEdgeBd; n<nEdge; n++) edge(n,5) = edge(n,4);
  
  nflag.set(-1);
  for (int n=nEdge-nEdgeBd; n<nEdge; n++){
    n1          = edge(n,0);
    n2          = edge(n,1);
    nflag(n1,0) = 0;
    nflag(n2,0) = 0;
  }

  for (int n=0; n<nEdge-nEdgeBd; n++){
    n1 = edge(n,0);
    n2 = edge(n,1);
    dx = x(n2,0)-x(n1,0);
    dy = x(n2,1)-x(n1,1);
    ds = 1./sqrt(dx*dx+dy*dy);
    dx = dx*ds;
    dy = dy*ds;
    if (nflag(n1,0) < 0)
      for (int m=0; m<csp2[n1]; m++){ //look at the cells surrounding the node
	j = csp1[n1][m];
	for (int k=0; k<3; k++) //find the 1st opposing node in the triangle
	  if (n1 != tri(j,k)){
	    k1 = tri(j,k);
	    break;
	  }
	for (int k=0; k<3; k++) //find the 2nd opposing node in the triangle
	  if (n1 != tri(j,k) && k1 != tri(j,k)){
	    k2 = tri(j,k);
	    break;
	  }
	dx1 = x(k1,0)-x(n1,0);
	dy1 = x(k1,1)-x(n1,1);
	ds  = 1./sqrt(dx1*dx1+dy1*dy1);
	dx1 = dx1*ds;
	dy1 = dy1*ds;
	dx2 = x(k2,0)-x(n1,0);
	dy2 = x(k2,1)-x(n1,1);
	ds  = 1./sqrt(dx2*dx2+dy2*dy2);
	dx2 = dx2*ds;
	dy2 = dy2*ds;
	cp1 = dx*dy1-dy*dx1;
	cp2 = dx*dy2-dy*dx2;
	if (cp1*cp2 <= 0. && j != edge(n,4) && j != edge(n,5)){
	  edge(n,4) = j;
	  break;
	}}
    if (nflag(n2,0) < 0)
      for (int m=0; m<csp2[n2]; m++){ //look at the cells surrounding the node
	j = csp1[n2][m];
	for (int k=0; k<3; k++) //find the 1st opposing node in the triangle
	  if (n2 != tri(j,k)){
	    k1 = tri(j,k);
	    break;
	  }
	for (int k=0; k<3; k++) //find the 2nd opposing node in the triangle
	  if (n2 != tri(j,k) && k1 != tri(j,k)){
	    k2 = tri(j,k);
	    break;
	  }
	dx1 = x(k1,0)-x(n2,0);
	dy1 = x(k1,1)-x(n2,1);
	ds  = 1./sqrt(dx1*dx1+dy1*dy1);
	dx1 = dx1*ds;
	dy1 = dy1*ds;
	dx2 = x(k2,0)-x(n2,0);
	dy2 = x(k2,1)-x(n2,1);
	ds  = 1./sqrt(dx2*dx2+dy2*dy2);
	dx2 = dx2*ds;
	dy2 = dy2*ds;
	cp1 = dx*dy1-dy*dx1;
	cp2 = dx*dy2-dy*dx2;
	if (cp1*cp2 <= 0. && j != edge(n,4) && j != edge(n,5)){
	  edge(n,5) = j;
	  break;
	}}
  }

  for (int n=nEdge-nEdgeBd; n<nEdge; n++){
    n1 = edge(n,0);
    n2 = edge(n,1);
    dx = x(n2,0)-x(n1,0);
    dy = x(n2,1)-x(n1,1);
    ds = 1./sqrt(dx*dx+dy*dy);
    dx = dx*ds;
    dy = dy*ds;
    for (int m=0; m<csp2[n1]; m++){ //look at the cells surrounding the node
      j = csp1[n1][m];
      if (j != edge(n,4) && j != edge(n,5)){
	for (int k=0; k<3; k++) //find the 1st opposing node in the triangle
	  if (n1 != tri(j,k)){
	    k1 = tri(j,k);
	    break;
	  }
	for (int k=0; k<3; k++) //find the 2nd opposing node in the triangle
	  if (n1 != tri(j,k) && k1 != tri(j,k)){
	    k2 = tri(j,k);
	    break;
	  }
	if (nflag(k1,0) == 0 || nflag(k2,0) == 0){
	  edge(n,4) = j;
	  break;
	}}}
    for (int m=0; m<csp2[n2]; m++){ //look at the cells surrounding the node
      j = csp1[n2][m];
      if (j != edge(n,4) && j != edge(n,5)){
	for (int k=0; k<3; k++) //find the 1st opposing node in the triangle
	  if (n2 != tri(j,k)){
	    k1 = tri(j,k);
	    break;
	  }
	for (int k=0; k<3; k++) //find the 2nd opposing node in the triangle
	  if (n2 != tri(j,k) && k1 != tri(j,k)){
	    k2 = tri(j,k);
	    break;
	  }
	if (nflag(k1,0) == 0 || nflag(k2,0) == 0){
	  edge(n,5) = j;
	  break;
	}}}
  }

  /*
  for (int n=0; n<nEdge; n++)
    cout << n << " "
	 << edge(n,0) << " " << edge(n,1) << " "
	 << edge(n,2) << " " << edge(n,3) << " "
	 << edge(n,4) << " " << edge(n,5) << endl;
  exit(0);
  */


  // deallocate work arrays
  for (int n=0; n<nNode; n++)
    if (csp1[n]) delete [] csp1[n];
  delete [] csp1;
  delete [] csp2;
  edgeT.deallocate();
  nflag.deallocate();
}
