#include "Strand2dFCSPTurbSA.h"


void Strand2dFCSPTurbSA::outputSolution(const int& npts,
				  ofstream& ffile,
				  const string& outputVar,
				  const double* q,
				  const double* qa,
				  const double* e,
				  const double* r)
{
  // freestream values
  int iq,iqa;
  double r0,p0,u0,v0,t0,s0,nu0;
  double* rValue;
  rValue = solution.getRefValues();
  p0     = rValue[0];
  u0     = rValue[1];
  v0     = rValue[2];
  t0     = rValue[3];
  nu0    = rValue[4];
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

 // nu_tilde
 else if (!strcmp(outputVar.c_str(),"nu_tilde")){
    for (int n=0; n<npts; n++){
      iq  = n*nq;
      iqa = n*nqa;
      ffile << qa[iqa+4] << endl; 
   }
 }

 // eddy_visc
 else if (!strcmp(outputVar.c_str(),"eddy_visc")){
   double rnu,mu,chi,chi3,fv1,mut;
   for (int n=0; n<npts; n++){
     iq   = n*nq;
     iqa  = n*nqa;
     rnu  = q [iq +4];
     mu   = qa[iqa+5];
     chi  = rnu/mu;
     chi3 = chi*chi*chi;
     fv1  = chi3/(chi3+cv1*cv1*cv1);
     mut  = rnu*fv1;
     if (rnu < 0.) mut = 0.;
     ffile << mut/mu << endl; 
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


  // R4
  else if (!strcmp(outputVar.c_str(),"R5")){
    for (int n=0; n<npts; n++){
      iq = n*nq;
      ffile << r[iq+4] << endl;
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


  // E5
  else if (!strcmp(outputVar.c_str(),"E5")){
    for (int n=0; n<npts; n++){
      iq = n*nq;
      ffile << e[iq+4] << endl;
    }
  }


  else{
    cout << "\n*** Output variable not recognized in outputSolution.C ***"
	 << endl;
    exit(0);
  }
}
