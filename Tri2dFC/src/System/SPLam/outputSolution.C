#include "Tri2dFCSPLam.h"


void Tri2dFCSPLam::outputSolution(const int& npts,
				  ofstream& ffile,
				  const string& outputVar,
				  const double* q,
				  const double* qa,
				  const double* e,
				  const double* r)
{
  // freestream values
  int iq,iqa;
  double r0,p0,u0,v0,t0,s0;
  double* rValue;
  rValue = solution.getRefValues();
  p0     = rValue[0];
  u0     = rValue[1];
  v0     = rValue[2];
  t0     = rValue[3];
  r0     = p0/(rGas*t0);
  s0     = p0/pow(r0,gamma);


  // pressure
  if      (!strcmp(outputVar.c_str(),"P")){
    for (int n=0; n<npts; n++){
      iqa = n*nqa;
      ffile << qa[iqa] << endl;
    }
  }


  // temperature
  else if (!strcmp(outputVar.c_str(),"T")){
    for (int n=0; n<npts; n++){
      iqa = n*nqa;
      ffile << qa[iqa+3] << endl;
    }
  }


  // density
  else if (!strcmp(outputVar.c_str(),"rho")){
    for (int n=0; n<npts; n++){
      iq = n*nq;
      ffile << q[iq] << endl;
    }
  }


  // entropy
  else if (!strcmp(outputVar.c_str(),"entropy")){
    for (int n=0; n<npts; n++){
      iq  = n*nq;
      iqa = n*nqa;
      ffile << fabs(qa[iqa]/pow(q[iq],gamma)-s0) << endl;
    }
  }


  // velocity
  else if (!strcmp(outputVar.c_str(),"U")){
    for (int n=0; n<npts; n++){
      iqa = n*nqa;
      ffile << qa[iqa+1] << "\t" << qa[iqa+2] << "\t" << 0. << endl;
    }
  }


  // R1
  else if (!strcmp(outputVar.c_str(),"R1")){
    for (int n=0; n<npts; n++){
      iq = n*nq;
      ffile << r[iq] << endl;
    }
  }


  // R2
  else if (!strcmp(outputVar.c_str(),"R2")){
    for (int n=0; n<npts; n++){
      iq = n*nq;
      ffile << r[iq+1] << endl;
    }
  }


  // R3
  else if (!strcmp(outputVar.c_str(),"R3")){
    for (int n=0; n<npts; n++){
      iq = n*nq;
      ffile << r[iq+2] << endl;
    }
  }


  // R4
  else if (!strcmp(outputVar.c_str(),"R4")){
    for (int n=0; n<npts; n++){
      iq = n*nq;
      ffile << r[iq+3] << endl;
    }
  }


  // E1
  else if (!strcmp(outputVar.c_str(),"E1")){
    for (int n=0; n<npts; n++){
      iq = n*nq;
      ffile << e[iq] << endl;
    }
  }


  // E2
  else if (!strcmp(outputVar.c_str(),"E2")){
    for (int n=0; n<npts; n++){
      iq = n*nq;
      ffile << e[iq+1] << endl;
    }
  }


  // E3
  else if (!strcmp(outputVar.c_str(),"E3")){
    for (int n=0; n<npts; n++){
      iq = n*nq;
      ffile << e[iq+2] << endl;
    }
  }


  // E4
  else if (!strcmp(outputVar.c_str(),"E4")){
    for (int n=0; n<npts; n++){
      iq = n*nq;
      ffile << e[iq+3] << endl;
    }
  }


  else{
    cout << "\n*** Output variable not recognized in outputSolution.C ***"
	 << endl;
    exit(0);
  }
}
