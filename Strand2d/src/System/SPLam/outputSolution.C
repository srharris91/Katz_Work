#include "StrandSPLam.h"
#include "StrandSPLamStepHeader.h"


void StrandSPLam::outputSolution(ofstream& ffile,
				 const int& npts,
				 const double* q,
				 const double* qa)
{
  // freestream values
  int iq,iqa;
  double r0,p0,u0,v0,t0,s0;
  double* rValue;
  rValue       = solution.getRefValues();
  p0           = rValue[0];
  u0           = rValue[1];
  v0           = rValue[2];
  t0           = rValue[3];
  r0           = p0/(rGas*t0);
  s0           = p0/pow(r0,gamma);

  // variables header
  ffile << "<PointData Scalars=\"P T rho entropy\" Vectors=\"U\">"
	<< endl


    // pressure
	<< "<DataArray type=\"Float32\" Name=\"P\" format=\"ascii\">"
	<< endl;
  for (int n=0; n<npts; n++){
    iqa = n*nqa;
    ffile << qa[iqa] << endl;
  }
  ffile << "</DataArray>" << endl


    // temperature
	<< "<DataArray type=\"Float32\" Name=\"T\" format=\"ascii\">"
	<< endl;
  for (int n=0; n<npts; n++){
    iqa = n*nqa;
    ffile << qa[iqa+3] << endl;
  }
  ffile << "</DataArray>" << endl


    // density
	<< "<DataArray type=\"Float32\" Name=\"rho\" format=\"ascii\">"
	<< endl;
  for (int n=0; n<npts; n++){
    iq = n*nq;
    ffile << q[iq] << endl;
  }
  ffile << "</DataArray>" << endl

    // entropy
	<< "<DataArray type=\"Float32\" Name=\"entropy\" format=\"ascii\">"
	<< endl;
  for (int n=0; n<npts; n++){
    iq  = n*nq;
    iqa = n*nqa;
    ffile << fabs(qa[iqa]/pow(q[iq],gamma)-s0) << endl;
  }
  ffile << "</DataArray>" << endl


    // velocity
	<< "<DataArray type=\"Float32\" Name=\"U\" "
	<< "NumberOfComponents=\"3\" format=\"ascii\">"
	<< endl;
  for (int n=0; n<npts; n++){
    iqa = n*nqa;
    ffile << qa[iqa+1] << "\t" << qa[iqa+2] << "\t" << 0. << endl;
  }
  ffile << "</DataArray>" << endl


    // end the data
	<< "</PointData>" << endl;
}
