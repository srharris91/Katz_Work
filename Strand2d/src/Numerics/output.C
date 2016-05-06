#include "StrandBlockSolver.h"


void StrandBlockSolver::output(const int& step)
{
  // output the step header information
  if (standAlone == 1) sys->outputInitialize(step);
  else{
    stringstream outp;
    string outps;
    outp << "mkdir -p output." << step;
    outps = outp.str();
    const char* coutps = outps.c_str();
    system(coutps);
    stringstream outp2;
    outp2 << "rm -f output." << step << "/*part*";
    outps = outp2.str();
    coutps = outps.c_str();
    system(coutps);
  }


  // form a general unstructured array of points and cells
  // from non-clipped cells
  int i,j;
  int mglevel=0;
  i = nPstr+1;
  j = nPstr+2;
  Array2D<int> pflag(i,nNodes);
  Array2D<int> cflag(j,nFaces);

  for (int n=0; n<nNodes; n++)
    for (int j=0; j<nPstr+1; j++) pflag(j,n) =-1;

  for (int n=0; n<nFaces; n++)
    for (int j=0; j<nPstr+2; j++) cflag(j,n) =-1;

  int n1,n2,ifc,ic,ip,jm,k;
  k = 0;
  for (int n=0; n<nFaces-nGfaces; n++){
    n1 = face(0,n);
    n2 = face(1,n);
  for (int j=1; j<fClip(n)+1; j++){
    jm           = j-1;
    cflag(j ,n ) = k;
    pflag(j ,n1) = 0;
    pflag(j ,n2) = 0;
    pflag(jm,n1) = 0;
    pflag(jm,n2) = 0;
    k++;
  }
  }
  int nCells = k;

  k = 0;
  for (int n=0; n<nNodes; n++)
  for (int j=0; j<nPstr+1; j++)
    if (pflag(j,n) == 0){
      pflag(j,n) = k;
    k++;
    }
  int nPoints = k;

  Array2D<double> xP(2,nPoints);
  for (int n=0; n<nNodes; n++)
  for (int j=0; j<nPstr+1; j++){
    k = pflag(j,n);
    if (k >= 0){
      xP(0,k) = x(0,j,n);
      xP(1,k) = x(1,j,n);
    }
  }

  Array2D<int> cell(4,nCells);
  for (int n=0; n<nFaces-nGfaces; n++){
    n1         = face(0,n);
    n2         = face(1,n);
    for (int j=1; j<fClip(n)+1; j++){
      jm        = j-1;
      k         = cflag(j,n);
      cell(0,k) = pflag(jm,n1);
      cell(1,k) = pflag(jm,n2);
      cell(2,k) = pflag(j ,n2);
      cell(3,k) = pflag(j ,n1);
    }}


  // if solution output, form nodal solution data and output to file
  if (iSolnFile != 0){
    nodalQFlag = 0;
    nodalQ(mglevel);
    nodalQa(mglevel);
    Array2D<double> qPP (nq ,nPoints);
    Array2D<double> qaPP(nqa,nPoints);
    for (int n=0; n<nNodes; n++)
    for (int j=0; j<nPstr+1; j++){
      k = pflag(j,n);
      if (k >= 0) for (int m=0; m<nq ; m++) qPP (m,k) = qp (m,j,n);
      if (k >= 0) for (int m=0; m<nqa; m++) qaPP(m,k) = qap(m,j,n);
    }

    ofstream ffile;
    ffile.open ("soln_part0.vtu"); //eventually replace 0 with id number
    ffile.setf(ios::scientific);
    ffile.precision(14);
    ffile << "<?xml version=\"1.0\"?>"
	  << endl
	  <<  "<VTKFile type=\"UnstructuredGrid\" version=\"0.1"
	  <<"\" byte_order=\"LittleEndian\">"
	  << endl
	  << "<UnstructuredGrid>"
	  << endl
	  << "<Piece NumberOfPoints=\"" << nPoints <<"\" NumberOfCells=\""
	  << nCells << "\">"
	  << endl;

    sys->outputSolution(ffile,
			nPoints,
			&qPP(0,0),
			&qaPP(0,0));

    ffile << "<Points>" << endl
	  << "<DataArray type=\"Float32\" NumberOfComponents"
	  << "=\"3\" format=\"ascii\">" << endl;
    for (int n=0; n<nPoints; n++)
      ffile << xP(0,n) << "\t" << xP(1,n) << "\t" << 0. << endl;
    ffile << "</DataArray>" << endl
	  << "</Points>" << endl
	  << "<Cells>" << endl
	  << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
	  << endl;
    for (int n=0; n<nCells; n++)
      ffile << cell(0,n) << "\t"
	    << cell(1,n) << "\t"
	    << cell(2,n) << "\t"
	    << cell(3,n) << endl;
    ffile << "</DataArray>" << endl
	  <<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
	  << endl;
    int k = 0;
    for (int n=0; n<nCells; n++){
      k += 4;
      ffile << k << endl;
    }
    ffile << "</DataArray>" << endl
	  << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">"
	  << endl;
    for (int n=0; n<nCells; n++) ffile << 9 << endl;
    ffile << "</DataArray>" << endl
	  << "</Cells>" << endl
	  << "</Piece>" << endl
	  << "</UnstructuredGrid>" << endl
	  << "</VTKFile>" << endl;
    ffile.close();
    stringstream ss;
    string sss;
    ss << "mv soln_part0.vtu output." << step;//eventually replace 0 w/id number
    sss = ss.str();
    const char* csss = sss.c_str();
    system(csss);
    qPP.deallocate();
    qaPP.deallocate();
  }


  // if residual output, form cell residual data and output to file
  if (iResdFile != 0){
    Array2D<double> rC(nq,nCells);
    for (int n=0; n<nFaces; n++)
    for (int j=0; j<nPstr+2; j++){
      k = cflag(j,n);
      if (k >= 0) for (int m=0; m<nq; m++) rC(m,k) = q(m,j,n)-q0(m,j,n);
    }

    ofstream ffile;
    ffile.open ("resd_part0.vtu"); //eventually replace 0 with id number
    ffile.setf(ios::scientific);
    ffile.precision(14);
    ffile << "<?xml version=\"1.0\"?>"
	  << endl
	  <<  "<VTKFile type=\"UnstructuredGrid\" version=\"0.1"
	  <<"\" byte_order=\"LittleEndian\">"
	  << endl
	  << "<UnstructuredGrid>"
	  << endl
	  << "<Piece NumberOfPoints=\"" << nPoints <<"\" NumberOfCells=\""
	  << nCells << "\">"
	  << endl;

    stringstream rn;
    rn << "<CellData Scalars=\"";
    for (k=0; k<nq; k++) rn << "r" << k+1 << " ";
    rn << "\">";
    ffile << rn.str() << endl;

    for (k=0; k<nq; k++){
      stringstream sn;
      sn << k+1;
      ffile << "<DataArray type=\"Float32\" Name=\"r" << sn.str()
	    << "\" format=\"ascii\">"
	    << endl;
      for (int n=0; n<nCells; n++) ffile << rC(k,n) << endl;
      ffile << "</DataArray>" << endl;
    }

    ffile << "</CellData>" << endl
	  << "<Points>" << endl
	  << "<DataArray type=\"Float32\" NumberOfComponents"
	  << "=\"3\" format=\"ascii\">" << endl;
    for (int n=0; n<nPoints; n++)
      ffile << xP(0,n) << "\t" << xP(1,n) << "\t" << 0. << endl;
    ffile << "</DataArray>" << endl
	  << "</Points>" << endl
	  << "<Cells>" << endl
	  << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
	  << endl;
    for (int n=0; n<nCells; n++)
      ffile << cell(0,n) << "\t"
	    << cell(1,n) << "\t"
	    << cell(2,n) << "\t"
	    << cell(3,n) << endl;
    ffile << "</DataArray>" << endl
	  <<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
	  << endl;
    int k = 0;
    for (int n=0; n<nCells; n++){
      k += 4;
      ffile << k << endl;
    }
    ffile << "</DataArray>" << endl
	  << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">"
	  << endl;
    for (int n=0; n<nCells; n++) ffile << 9 << endl;
    ffile << "</DataArray>" << endl
	  << "</Cells>" << endl
	  << "</Piece>" << endl
	  << "</UnstructuredGrid>" << endl
	  << "</VTKFile>" << endl;
    ffile.close();
    stringstream ss;
    string sss;
    ss << "mv resd_part0.vtu output." << step;//eventually replace 0 w/id number
    sss = ss.str();
    const char* csss = sss.c_str();
    system(csss);
    rC.deallocate();
  }


  // if error output, form cell error data and output to file
  if (iErrFile != 0){
    Array2D<double> eC(nq,nCells);
    int npts =(nPstr+2)*(nFaces+nBedges);
    j = nPstr+2;
    k = nFaces+nBedges;
    Array3D<double> qe(nq,j,k);
    sys->initFlow(npts,
		  &xc(0,0,0),
		  &qe(0,0,0));

    int kk;
    for (int n=0; n<nFaces; n++)
    for (int j=0; j<nPstr+2; j++){
      k = cflag(j,n);
      if (k >= 0) for (int m=0; m<nq; m++) eC(m,k) = q(m,j,n)-qe(m,j,n);
    }

    ofstream ffile;
    ffile.open ("err_part0.vtu"); //eventually replace 0 with id number
    ffile.setf(ios::scientific);
    ffile.precision(14);
    ffile << "<?xml version=\"1.0\"?>"
	  << endl
	  <<  "<VTKFile type=\"UnstructuredGrid\" version=\"0.1"
	  <<"\" byte_order=\"LittleEndian\">"
	  << endl
	  << "<UnstructuredGrid>"
	  << endl
	  << "<Piece NumberOfPoints=\"" << nPoints <<"\" NumberOfCells=\""
	  << nCells << "\">"
	  << endl;

    stringstream rn;
    rn << "<CellData Scalars=\"";
    for (k=0; k<nq; k++) rn << "e" << k+1 << " ";
    rn << "\">";
    ffile << rn.str() << endl;

    for (k=0; k<nq; k++){
      stringstream sn;
      sn << k+1;
      ffile << "<DataArray type=\"Float32\" Name=\"e" << sn.str()
	    << "\" format=\"ascii\">"
	    << endl;
      for (int n=0; n<nCells; n++) ffile << eC(k,n) << endl;
      ffile << "</DataArray>" << endl;
    }

    ffile << "</CellData>" << endl
	  << "<Points>" << endl
	  << "<DataArray type=\"Float32\" NumberOfComponents"
	  << "=\"3\" format=\"ascii\">" << endl;
    for (int n=0; n<nPoints; n++)
      ffile << xP(0,n) << "\t" << xP(1,n) << "\t" << 0. << endl;
    ffile << "</DataArray>" << endl
	  << "</Points>" << endl
	  << "<Cells>" << endl
	  << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
	  << endl;
    for (int n=0; n<nCells; n++)
      ffile << cell(0,n) << "\t"
	    << cell(1,n) << "\t"
	    << cell(2,n) << "\t"
	    << cell(3,n) << endl;
    ffile << "</DataArray>" << endl
	  <<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
	  << endl;
    int k = 0;
    for (int n=0; n<nCells; n++){
      k += 4;
      ffile << k << endl;
    }
    ffile << "</DataArray>" << endl
	  << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">"
	  << endl;
    for (int n=0; n<nCells; n++) ffile << 9 << endl;
    ffile << "</DataArray>" << endl
	  << "</Cells>" << endl
	  << "</Piece>" << endl
	  << "</UnstructuredGrid>" << endl
	  << "</VTKFile>" << endl;
    ffile.close();
    stringstream ss;
    string sss;
    ss << "mv err_part0.vtu output." << step;//eventually replace 0 w/id number
    sss = ss.str();
    const char* csss = sss.c_str();
    system(csss);

    ofstream efile;
    efile.open("error.dat",ios::app);
    efile.setf(ios::scientific);
    efile.precision(14);
    double erms[nq];
    for (int k=0; k<nq; k++) erms[k] = 0.;
    for (int n=0; n<nCells; n++)
      for (int k=0; k<nq; k++)
	erms[k] +=(eC(k,n)*eC(k,n));
    for (int k=0; k<nq; k++) erms[k] = sqrt(erms[k]/((double)nCells));
    efile << sqrt((double)nCells) << "\t";
    for (int k=0; k<nq; k++) efile << erms[k] << "\t";
    efile << endl;
    efile.close();
    cout << "\nError statistics: " << endl;
    cout << sqrt((double)nCells) << "\t";
    for (int k=0; k<nq; k++) cout << erms[k] << "\t";
    cout << endl;

    eC.deallocate();
    qe.deallocate();
  }


  // if surface data, form surface data and output to file
  if (iSurfFile != 0){
    nodalQFlag = 0;
    nodalQ(mglevel);
    ofstream ffile;
    ffile.open ("surf_data.dat");
    ffile.setf(ios::scientific);
    ffile.precision(14);
    int j = 0;
    int ixc,iq;
    double xn,yn;
    for (int n=0; n<nNodes-nGnodes; n++){
      xn = x(0,j,n);
      yn = x(1,j,n);
      ffile << xn << "\t" << yn << "\t";
      for (int k=0; k<nq; k++) ffile << qp(k,j,n) << "\t";
      ffile << endl;
    }
    ffile.close();

    // output surface forces
    forcex = 0.;
    forcey = 0.;
    j = 0;
    int npts = 1;
    int ifu,iqa,iqx,iqax;
    Array2D<double> force(ndim,npts);
    for (int n=0; n<nFaces-nGfaces; n++){
      sys->outputSurfaceForces(npts,
			       &fTag(n),
			       &facu(0,j,n),
			       &q(0,j,n),
			       &qa(0,j,n),
			       &qx(0,0,j,n),
			       &qax(0,0,j,n),
			       &force(0,0));
      forcex += force(0,0);
      forcey += force(1,0);
    }
    cout << "\nx-component force: " << forcex << endl;
    cout <<   "y-component force: " << forcey << endl;
    force.deallocate();
  }


  // deallocate variables
  pflag.deallocate();
  cflag.deallocate();
  xP.deallocate();
  cell.deallocate();
}
