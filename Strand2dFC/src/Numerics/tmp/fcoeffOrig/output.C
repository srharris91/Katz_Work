#include "Strand2dFCBlockSolver.h"
#include "lagrangePoly.h"


void Strand2dFCBlockSolver::output(const int& nBlocks,
				const int& step)
{
  // create a fully unstructured grid from the strand grid
  int nQuad=nSurfElem*meshOrder*(nStrandNode-1);
  int n1,n2;
  double x0,y0,nx,ny;
  Array3D<double> x(nSurfNode,nStrandNode,2);
  Array2D<int> map(nSurfNode,nStrandNode);
  Array2D<int> quad(nQuad,4);
  Array2D<int> surfCell(meshOrder,2);

  if (meshOrder == 1){
    surfCell(0,0) = 0;
    surfCell(0,1) = 1;
  }
  else if (meshOrder == 2){
    surfCell(0,0) = 0;
    surfCell(0,1) = 2;
    surfCell(1,0) = 2;
    surfCell(1,1) = 1;
  }
  else if (meshOrder == 3){
    surfCell(0,0) = 0;
    surfCell(0,1) = 2;
    surfCell(1,0) = 2;
    surfCell(1,1) = 3;
    surfCell(2,0) = 3;
    surfCell(2,1) = 1;
  }

  int k=0;
  for (int n=0; n<nSurfNode; n++){
    x0 = surfX(n,0);
    y0 = surfX(n,1);
    nx = pointingVec(n,0);
    ny = pointingVec(n,1);
    for (int j=0; j<nStrandNode; j++){
      x(n,j,0) = x0+nx*strandX(j);
      x(n,j,1) = y0+ny*strandX(j);
      map(n,j) = k++;
    }}

  nQuad = 0;
  for (int n=0; n<nSurfElem; n++)
    for (int i=0; i<meshOrder; i++){
      n1 = surfElem(n,surfCell(i,0));
      n2 = surfElem(n,surfCell(i,1));
      for (int j=0; j<nStrandNode-1; j++){
	quad(nQuad  ,0) = map(n1,j  );
	quad(nQuad  ,1) = map(n2,j  );
	quad(nQuad  ,2) = map(n2,j+1);
	quad(nQuad++,3) = map(n1,j+1);
      }}

  // create output directory for this unsteady step
  stringstream a;
  a.str("");
  a.clear();
  a << "mkdir -p output." << step;
  system(a.str().c_str());
  a.str("");
  a.clear();


  // output solution file
  if (nOutputVars > 0){

    // determine error at the plotting points
    Array3D<double> qe(nSurfNode,nStrandNode,nq),er(nSurfNode,nStrandNode,nq);
    sys->initFlow(nSurfNode*nStrandNode,
                  &x(0,0,0),
                  &qe(0,0,0));

    for (int n=0; n<nSurfNode; n++)
      for (int j=0; j<nStrandNode; j++)
	for (int k=0; k<nq; k++) er(n,j,k) =(q(n,j,k)-qe(n,j,k))/rmsNorm(k);


    // output the step header information if on block 0
    if (ID == 0){
      ofstream ffile;
      a << "strandSolution" << step << ".pvtu";
      ffile.open (a.str().c_str());
      a.str("");
      a.clear();
      ffile.setf(ios::scientific);
      ffile.precision(14);

      ffile << "<?xml version=\"1.0\"?>"
            << endl
            << "<VTKFile type=\"PUnstructuredGrid\" "
            << "version=\"0.1\" byte_order=\"LittleEndian\">"
            << endl
            << "<PUnstructuredGrid GhostLevel=\"0\">"
            << endl
            << "<PPointData ";

      if (nOutputScalars > 0){
        ffile << "Scalars=\"";
        int k=0;
        for (int n=0; n<nOutputVars; n++)
          if (outputVarLength(n) == 1){
            if (k != 0) ffile << " ";
            k++;
            ffile << outputVars(n);
          }
        ffile << "\"";
        if (nOutputVectors > 0) ffile << " ";
      }

      if (nOutputVectors > 0){
        ffile << "Vectors=\"";
        int k=0;
        for (int n=0; n<nOutputVars; n++)
          if (outputVarLength(n) == 3){
            if (k != 0) ffile << " ";
            k++;
            ffile << outputVars(n);
          }
        ffile << "\"";
      }

      ffile << ">"
            << endl;

      for (int n=0; n<nOutputVars; n++){
        ffile << "<PDataArray type=\"Float32\" Name=\""
              << outputVars(n) << "\"";
        if (outputVarLength(n) == 3) ffile << " NumberOfComponents=\"3\"";
        ffile << "/>"
              << endl;
      }

      ffile << "</PPointData>"
            << endl
            << "<PPoints>"
            << endl
            << "<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>"
            << endl
            << "</PPoints>"
            << endl;

      for (int n=0; n<nBlocks; n++){
        a << "<Piece Source=\"output." << step
          << "/strandSolutionBlock" << n << ".vtu\"/>";
        ffile << a.str()
              << endl;
        a.str("");
        a.clear();
      }

      ffile << "</PUnstructuredGrid>"
            << endl
            << "</VTKFile>"
            << endl;

      ffile.close();
    }


    // write data to file for this block
    ofstream ffile;
    a << "strandSolutionBlock" << ID << ".vtu";
    ffile.open (a.str().c_str());
    ffile.setf(ios::scientific);
    ffile.precision(14);
    a.str("");
    a.clear();

    ffile << "<?xml version=\"1.0\"?>"
          << endl
          <<  "<VTKFile type=\"UnstructuredGrid\" version=\"0.1"
          <<"\" byte_order=\"LittleEndian\">"
          << endl
          << "<UnstructuredGrid>"
          << endl
          << "<Piece NumberOfPoints=\""
          << nSurfNode*nStrandNode
          <<"\" NumberOfCells=\""
          << nQuad << "\">"
          << endl
          << "<PointData ";
    if (nOutputScalars > 0){
      ffile << "Scalars=\"";
      int k=0;
      for (int n=0; n<nOutputVars; n++)
        if (outputVarLength(n) == 1){
          if (k != 0) ffile << " ";
          k++;
          ffile << outputVars(n);
        }
      ffile << "\"";
      if (nOutputVectors > 0) ffile << " ";
    }

    if (nOutputVectors > 0){
      ffile << "Vectors=\"";
      int k=0;
      for (int n=0; n<nOutputVars; n++)
        if (outputVarLength(n) == 3){
          if (k != 0) ffile << " ";
          k++;
          ffile << outputVars(n);
        }
      ffile << "\"";
    }

    ffile << ">"
          << endl;

    for (int n=0; n<nOutputVars; n++){
      if (outputVarLength(n) == 1){
        ffile << "<DataArray type=\"Float32\" Name=\""
              << outputVars(n)
              << "\" format=\"ascii\">"
              << endl;
      }
      else{
        ffile << "<DataArray type=\"Float32\" Name=\""
              << outputVars(n)
              << "\" NumberOfComponents=\"3\" format=\"ascii\">"
              << endl;
      }
      // get rid of very small numbers the Paraview doesn't like
      double eps=1.e-14;
      for (int nn=0; nn<nSurfNode; nn++)
	for (int j=0; j<nStrandNode; j++){
	  for (int k=0; k<nq; k++)
	    if (fabs(q(nn,j,k)) < eps) q(nn,j,k) = 0.;
	  for (int k=0; k<nqa; k++)
	    if (fabs(qa(nn,j,k)) < eps) qa(nn,j,k) = 0.;
	  for (int k=0; k<nq; k++)
	    if (fabs(er(nn,j,k)) < eps) er(nn,j,k) = 0.;
	  for (int k=0; k<nq; k++)
	    if (fabs(r(nn,j,k)) < eps) r(nn,j,k) = 0.;
	}

      sys->outputSolution(nSurfNode*nStrandNode,
                          ffile,
                          outputVars(n),
                          &q(0,0,0),
                          &qa(0,0,0),
                          &er(0,0,0),
                          &r(0,0,0));
      ffile << "</DataArray>" << endl;
    }

    ffile << "</PointData>" << endl
          << "<Points>" << endl
          << "<DataArray type=\"Float32\" NumberOfComponents"
          << "=\"3\" format=\"ascii\">" << endl;
    double xn,yn,eps=1.e-14;
    for (int n=0; n<nSurfNode; n++)
      for (int j=0; j<nStrandNode; j++){
	xn = x(n,j,0);
	yn = x(n,j,1);
	if (fabs(xn) < eps) xn = 0.;
	if (fabs(yn) < eps) yn = 0.;
	ffile << xn << "\t" << yn << "\t" << 0. << endl;
      }
    ffile << "</DataArray>" << endl
          << "</Points>" << endl
          << "<Cells>" << endl
          << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
          << endl;
    for (int n=0; n<nQuad; n++)
      ffile << quad(n,0) << "\t"
            << quad(n,1) << "\t"
            << quad(n,2) << "\t"
            << quad(n,3) << endl;
    ffile << "</DataArray>" << endl
          <<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
          << endl;
    int k = 0;
    for (int n=0; n<nQuad; n++){
      k += 4;
      ffile << k << endl;
    }
    ffile << "</DataArray>" << endl
          << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">"
          << endl;
    for (int n=0; n<nQuad; n++) ffile << 9 << endl;
    ffile << "</DataArray>" << endl
          << "</Cells>" << endl
          << "</Piece>" << endl
          << "</UnstructuredGrid>" << endl
          << "</VTKFile>" << endl;
    ffile.close();

    a << "mv strandSolutionBlock" << ID << ".vtu output." << step;
    system(a.str().c_str());
    a.str("");
    a.clear();


    // output L2 norm of solution error
    if (iErrFile != 0){
      ofstream efile;
      //a << "strandErrorBlock" << ID << ".dat";
      a << "error.dat";
      efile.open(a.str().c_str(),ios::app);
      efile.setf(ios::scientific);
      efile.precision(14);
      a.str("");
      a.clear();

      double erms[nq];
      for (int k=0; k<nq; k++) erms[k] = 0.;
      for (int n=0; n<nSurfNode; n++)
	for (int j=0; j<nStrandNode; j++)
	  for (int k=0; k<nq; k++) erms[k] += pow(er(n,j,k),2);
      for (int k=0; k<nq; k++)
	erms[k] = sqrt(erms[k]/(double)(nSurfNode*nStrandNode));
      efile << sqrt((double)(nSurfNode*nStrandNode)) << "\t";
      for (int k=0; k<nq; k++) efile << erms[k] << "\t";
      efile << endl;
      efile.close();

      //a << "mv strandErrorBlock" << ID << ".dat output." << step;
      //system(a.str().c_str());
      //a.str("");
      //a.clear();

      cout.setf(ios::scientific);
      //cout << "\nError statistics: " << endl;
      //cout << sqrt((double)(nSurfNode*nStrandNode)) << "\t";
      //for (int k=0; k<nq; k++) cout << erms[k] << "\t";
      //cout << endl;
    }

    // deallocate work arrays
    qe.deallocate();
    er.deallocate();
  }


  // output surface data to file
  if (iSurfFile != 0){
    ofstream sfile;
    a << "forces.dat";
    if (step == 0) sfile.open (a.str().c_str());
    else sfile.open (a.str().c_str(),ios::app);
    sfile.setf(ios::scientific);
    sfile.precision(14);
    a.str("");
    a.clear();

    // obtain solution gradients at the nodes
    int tag,ni,nm;
    double jQ,sxQ,nxQ,syQ,nyQ,qsQ[nq],qasQ[nqa],qnQ[nq],qanQ[nqa],
      qQ[nq],qaQ[nqa],qxQ[nq],qyQ[nq],qaxQ[nqa],qayQ[nqa],force[ndim],
      forcex=0.,forcey=0.,dnr=1./deltaN;
    Array2D<double> qnN(meshOrder+1,nq),qanN(meshOrder+1,nqa);
    for (int n=0; n<nSurfElem; n++){
      tag = surfElemTag(n); // boundary tag for the element

      // gradients in the n-direction at nodes in this element
      qnN.set(0.);
      qanN.set(0.);
      for (int i=0; i<meshOrder+1; i++){
	ni = surfElem(n,i);
	for (int m=0; m<nIcb; m++){
	  for (int k=0; k<nq ; k++) qnN (i,k) += dnr*icb(0,m)*q (ni,m,k);
	  for (int k=0; k<nqa; k++) qanN(i,k) += dnr*icb(0,m)*qa(ni,m,k);
	}}

      // add the contributions from each quadrature point in this element
      for (int i=0; i<nQuadPoint; i++){

	// gradients in the s-direction at quadrature point
	for (int k=0; k<nq ; k++) qsQ [k] = 0.;
	for (int k=0; k<nqa; k++) qasQ[k] = 0.;
        for (int m=0; m<meshOrder+1; m++){ // mth Lagrange poly. in mapping
          nm = surfElem(n,m);
	  for (int k=0; k<nq ; k++) qsQ [k] += lsQ(i,m)*q (nm,0,k);
	  for (int k=0; k<nqa; k++) qasQ[k] += lsQ(i,m)*qa(nm,0,k);
        }

	// interpolate n-direction gradients to quadrature point
	for (int k=0; k<nq ; k++) qnQ [k] = 0.;
	for (int k=0; k<nqa; k++) qanQ[k] = 0.;
        for (int m=0; m<meshOrder+1; m++){ // mth Lagrange poly. in mapping
	  for (int k=0; k<nq ; k++) qnQ [k] += lQ(i,m)*qnN (m,k);
	  for (int k=0; k<nqa; k++) qanQ[k] += lQ(i,m)*qanN(m,k);
        }

	// x-y gradients at quadrature point
	jQ  = 1./jacQ(n,i);
	sxQ = ynQ(n,i)*jQ;
	nxQ =-ysQ(n,i)*jQ;
	syQ =-xnQ(n,i)*jQ;
	nyQ = xsQ(n,i)*jQ;
	for (int k=0; k<nq ; k++) qxQ [k] = qsQ [k]*sxQ+qnQ [k]*nxQ;
	for (int k=0; k<nq ; k++) qyQ [k] = qsQ [k]*syQ+qnQ [k]*nyQ;
	for (int k=0; k<nqa; k++) qaxQ[k] = qasQ[k]*sxQ+qanQ[k]*nxQ;
	for (int k=0; k<nqa; k++) qayQ[k] = qasQ[k]*syQ+qanQ[k]*nyQ;

	// interpolate q to quadrature point
	for (int k=0; k<nq ; k++) qQ [k] = 0.;
	for (int k=0; k<nqa; k++) qaQ[k] = 0.;
        for (int m=0; m<meshOrder+1; m++){ // mth Lagrange poly. in mapping
          nm = surfElem(n,m);
	  for (int k=0; k<nq ; k++) qQ [k] += lQ(i,m)*q (nm,0,k);
	  for (int k=0; k<nqa; k++) qaQ[k] += lQ(i,m)*qa(nm,0,k);
        }

	// compute surface force at quadradure point
        sys->outputSurfaceForces(1,
                                 &tag,
                                 &xsQ(n,i),
				 &ysQ(n,i),
                                 &qQ[0],
                                 &qaQ[0],
                                 &qxQ[0],
                                 &qyQ[0],
                                 &qaxQ[0],
                                 &qayQ[0],
                                 &force[0]);

	// add the surface force contribution with quadrature weight
        forcex += wQ(i)*force[0];
        forcey += wQ(i)*force[1];
      }}

    // output surface forces to screen and file
    cout << "\nx-component force: " << forcex << endl;
    cout <<   "y-component force: " << forcey << endl;
    cout << "\n";
    sfile << step << " " << double(step)*dtUnsteady
          << " " << forcex << " " << forcey << endl;
    sfile.close();
    qnN.deallocate();
    qanN.deallocate();
  }


  // clean up
  x.deallocate();
  map.deallocate();
  quad.deallocate();
  surfCell.deallocate();
}
