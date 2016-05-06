#include "Perturb.h"

///\cond extern
extern "C"
{
  void inputsolutionreadperturb_(const int&,
				 const char*,
				 const int&,
				 double&,
				 double*);
}
///\endcond


// [Perturb]
Perturb::Perturb()
{
}
// [Perturb]


// [~Perturb]
Perturb::~Perturb()
{
}
// [~Perturb]


// [initialize]
void Perturb::initialize(const int& nq0,
			 const int& ndim0,
			 const string& inputFile,
			 State* state0,
			 Transport* transport0)
{
  nq = nq0;
  ndim = ndim0;
  rValue = new double[nq];
  inputsolutionreadperturb_(inputFile.size(),
			    inputFile.c_str(),
			    nq,
			    pert,
			    rValue);
  state     = state0;
  transport = transport0;
}
void Perturb::initialize(const int& nq0,
			 const int& ndim0,
			 const string& inputFile)
{
  nq = nq0;
  ndim = ndim0;
  rValue = new double[nq];
  inputsolutionreadperturb_(inputFile.size(),
			    inputFile.c_str(),
			    nq,
			    pert,
			    rValue);
}
// [initialize]


// [getQ]
void Perturb::getQ(const int& npts,
		   const double* x,
		   double* q)
{
  int i;
  double a;
  for (int n=0; n<npts; n++){
    i      = n*nq;
  for (int k=0; k<nq; k++){
    a      = double(rand())/double(RAND_MAX);
    a      =(a-.5)*2.*pert; // random number between -.1 and .1
    q[i+k] = rValue[k]*(1.+a);
  }}
}
// [getQ]


// [finalize]
void Perturb::finalize()
{
}
// [finalize]
