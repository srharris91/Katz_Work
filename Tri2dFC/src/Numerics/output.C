#include "Tri2dFCBlockSolver.h"
#include "lagrangePoly.h"


void Tri2dFCBlockSolver::output(const int& nBlocks,
				const int& step)
{
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
    Array2D<double> qe(nNode,nq),er(nNode,nq);
    sys->initFlow(nNode,
		  &x(0,0),
		  &qe(0,0));

    for (int n=0; n<nNode; n++)
      for (int k=0; k<nq; k++) er(n,k) =(q(n,k)-qe(n,k))/rmsNorm(k);


    // output the step header information if on block 0
    if (ID == 0){
      ofstream ffile;
      a << "triSolution" << step << ".pvtu";
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
	  << "/triSolutionBlock" << n << ".vtu\"/>";
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
    a << "triSolutionBlock" << ID << ".vtu";
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
	  << nNode
	  <<"\" NumberOfCells=\""
	  << nTri << "\">"
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
      sys->outputSolution(nNode,
			  ffile,
			  outputVars(n),
			  &q(0,0),
			  &qa(0,0),
			  &er(0,0),
			  &r(0,0));
      ffile << "</DataArray>" << endl;
    }
    
    ffile << "</PointData>" << endl
	  << "<Points>" << endl
	  << "<DataArray type=\"Float32\" NumberOfComponents"
	  << "=\"3\" format=\"ascii\">" << endl;
    double xn,yn,eps=1.e-14;
    for (int n=0; n<nNode; n++){
      xn = x(n,0);
      yn = x(n,1);
      if (fabs(x(n,0)) < eps) xn = 0.;
      if (fabs(x(n,1)) < eps) yn = 0.;
      ffile << xn << "\t" << yn << "\t" << 0. << endl;
    }
    ffile << "</DataArray>" << endl
	  << "</Points>" << endl
	  << "<Cells>" << endl
	  << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
	  << endl;
    for (int n=0; n<nTri; n++)
      ffile << tri(n,0) << "\t"
	    << tri(n,1) << "\t"
	    << tri(n,2) << endl;
    ffile << "</DataArray>" << endl
	  <<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
	  << endl;
    int k = 0;
    for (int n=0; n<nTri; n++){
      k += 3;
      ffile << k << endl;
    }
    ffile << "</DataArray>" << endl
	  << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">"
	  << endl;
    for (int n=0; n<nTri; n++) ffile << 5 << endl;
    ffile << "</DataArray>" << endl
	  << "</Cells>" << endl
	  << "</Piece>" << endl
	  << "</UnstructuredGrid>" << endl
	  << "</VTKFile>" << endl;
    ffile.close();

    a << "mv triSolutionBlock" << ID << ".vtu output." << step;
    system(a.str().c_str());
    a.str("");
    a.clear();


    // output L2 norm of solution error
    if (iErrFile != 0){
      ofstream efile;
      //a << "triErrorBlock" << ID << ".dat";
      a << "error.dat";
      efile.open(a.str().c_str(),ios::app);
      efile.setf(ios::scientific);
      efile.precision(14);
      a.str("");
      a.clear();

      double erms[nq];
      for (int k=0; k<nq; k++) erms[k] = 0.;
      for (int n=0; n<nNode; n++)
	for (int k=0; k<nq; k++) erms[k] += pow(er(n,k),2);
      for (int k=0; k<nq; k++) erms[k] = sqrt(erms[k]/(double)nNode);
      efile << sqrt((double)nNode) << "\t";
      for (int k=0; k<nq; k++) efile << erms[k] << "\t";
      efile << endl;
      efile.close();

      //a << "mv triErrorBlock" << ID << ".dat output." << step;
      //system(a.str().c_str());
      //a.str("");
      //a.clear();

      cout.setf(ios::scientific);
      cout << "\nError statistics: " << endl;
      cout << sqrt((double)nNode) << "\t";
      for (int k=0; k<nq; k++) cout << erms[k] << "\t";
      cout << endl;
    }


    // deallocate work arrays
    qe.deallocate();
    er.deallocate();
  }


  // output surface data to file
  if (iSurfFile != 0){
    ofstream sfile;
    a << "triSurfaceBlock" << ID << ".dat";
    if (step == 0) sfile.open (a.str().c_str());
    else sfile.open (a.str().c_str(),ios::app);
    sfile.setf(ios::scientific);
    sfile.precision(14);
    a.str("");
    a.clear();

    Array3D<double> qax;
    qax.allocate(nNode,2,nqa);
    gradient(nq,&q(0,0),&qx(0,0,0));
    gradient(nqa,&qa(0,0),&qax(0,0,0));

    int m=0;
    for (int n=nNode-nNodeBd; n<nNode; n++)
      sys->outputSurfaceSolution(1,
				 sfile,
				 &nodeBd(m++),
				 &x(n,0),
				 &x(n,1),
				 &q(n,0),
				 &qa(n,0),
				 &qx(n,0,0),
				 &qx(n,1,0),
				 &qax(n,0,0),
				 &qax(n,1,0));

    qax.deallocate();

    a << "mv triSurfaceBlock" << ID << ".dat output." << step;
    system(a.str().c_str());
    a.str("");
    a.clear();

    ofstream ffile; //force file
    a << "forces.dat";
    if (step == 0) ffile.open (a.str().c_str());
    else ffile.open (a.str().c_str(),ios::app);
    ffile.setf(ios::scientific);
    ffile.precision(14);
    a.str("");
    a.clear();

    forcex = 0.;
    forcey = 0.;
    int jj,nn,tag;
    double dd,dx,dy,
      qQ[nq],qxQ[nq],qyQ[nq],qaQ[nqa],qaxQ[nqa],qayQ[nqa],force[ndim];
    for (int n=0; n<nEdgeQ; n++){
      nn  = edgeQ(n,0); // element on which lies this edge
      jj  = edgeQ(n,1); // edge number on the element
      tag = edgeQ(n,2); // boundary tag for the edge
      for (int i=0; i<nsq; i++){ //for each surface quadrature point

	// interpolate solution and gradients to the quadrature points
	for (int k=0; k<nq ; k++) qQ[k]   = 0.;
	for (int k=0; k<nq ; k++) qxQ[k]  = 0.;
	for (int k=0; k<nq ; k++) qyQ[k]  = 0.;
	for (int k=0; k<nqa; k++) qaQ[k]  = 0.;
	for (int k=0; k<nqa; k++) qaxQ[k] = 0.;
	for (int k=0; k<nqa; k++) qayQ[k] = 0.;
	for (int j=0; j<nne; j++){
	  dd = gxS(n,i,j,0);
	  dx = gxS(n,i,j,1);
	  dy = gxS(n,i,j,2);
	  m  = elem(nn,j);
	  for (int k=0; k<nq ; k++){
	    qQ  [k] += q (m,k)*dd;
	    qxQ [k] += q (m,k)*dx;
	    qyQ [k] += q (m,k)*dy;
	  }
	  for (int k=0; k<nqa; k++){
	    qaQ [k] += qa(m,k)*dd;
	    qaxQ[k] += qa(m,k)*dx;
	    qayQ[k] += qa(m,k)*dy;
	  }}

	// compute surface force contribution
	sys->outputSurfaceForces(1,
				 &tag,
				 &nxQ(n,i,0),
				 &qQ[0],
				 &qxQ[0],
				 &qyQ[0],
				 &qaQ[0],
				 &qaxQ[0],
				 &qayQ[0],
				 &force[0]);
	forcex += force[0];
	forcey += force[1];
      }}

    cout << "\nx-component force: " << forcex << endl;
    cout <<   "y-component force: " << forcey << endl;
    cout << "\n";
    ffile << step << " " << double(step)*dtUnsteady
	  << " " << forcex << " " << forcey << endl;
    ffile.close();
  }
}
