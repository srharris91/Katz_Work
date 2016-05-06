#include "Strand2dFCSPLamBc.h"


// [Strand2dFCSPLamBc]
Strand2dFCSPLamBc::Strand2dFCSPLamBc()
{
  nq = 0;
  nqa = 0;
  ncomp = 0;
  ndim = 0;
  gamma = 0.;
  gm1 = 0.;
  ggm1 = 0.;
  rGas = 0.;
  bValue = NULL;
  state = NULL;
  transport = NULL;
  solution = NULL;
}
// [Strand2dFCSPLamBc]


// [~Strand2dFCSPLamBc]
Strand2dFCSPLamBc::~Strand2dFCSPLamBc()
{
}
// [~Strand2dFCSPLamBc]


// [initialize]
void Strand2dFCSPLamBc::initialize(const int& nq0,
				const int& nqa0,
				const int& ndim0,
				const int& ncomp0,
				const int& viscous0,
				const double& gamma0,
				const double& rGas0,
				double* bValue0,
				State* state0,
				Transport* transport0,
				Solution* solution0)
{
  nq = nq0;
  nqa = nqa0;
  ndim = ndim0;
  ncomp = ncomp0;
  viscous = viscous0;
  gamma = gamma0;
  rGas = rGas0;
  gm1 = gamma-1.;
  ggm1 = gamma/gm1;
  if (!bValue) bValue = new double[nq];
  for (int n=0; n<nq; n++) bValue[n] = bValue0[n];
  state = state0;
  transport = transport0;
  solution = solution0;
  transport->getPrn(&Prn);
}
// [initialize]


// [BCVector]
void Strand2dFCSPLamBc::BCVector(const double* nx,
				 const double* q,
				 const double* qa,
				 double* rb)
{
  for (int n=0; n<nq; n++) rb[n] = 0.; //none
}
// [BCVector]


// [BCVectorSelfJacobian]
void Strand2dFCSPLamBc::BCVectorSelfJacobian(const double* nx,
					     const double* q,
					     const double* qa,
					     double* M)
{
  for (int n=0; n<nq*nq; n++) M[n] = 0.;
}
// [BCVectorSelfJacobian]


// [BCSelectionMatrix]
void Strand2dFCSPLamBc::BCSelectionMatrix(const double* nx,
					  const double* q,
					  const double* qa,
					  double* L)
{
  for (int n=0; n<nq*nq; n++) L[n     ] = 0.;
  for (int n=0; n<nq   ; n++) L[n*nq+n] = 1.; //all EOMs
}
// [BCSelectionMatrix]


// [BCPenalty]
void Strand2dFCSPLamBc::BCPenalty(const int& inout,
				  const double* A,
				  const double& Pinv0,
				  const double* q,
				  const double* qa,
				  const double* g,
				  const double* uw,
				  double* rb)
{
  for (int n=0; n<nq; n++) rb[n] = 0.; 
}
// [BCPenalty]


// [BCPenaltyVis]
void Strand2dFCSPLamBc::BCPenaltyVis(const double& Pinv0,
				     const double* q,
				     const double* qa,
				     const double* gv,
				     const double* uw,
				     double* rb)
{
  for (int n=0; n<nq; n++) rb[n] = 0.; 
}
// [BCPenaltyVis]


// [BCPenaltyJacobian]
void Strand2dFCSPLamBc::BCPenaltyJacobian(const int& inout,
					  const double* A,
					  const double& Pinv0,
					  const double* q,
					  const double* qa,
					  const double* g,
					  const double* uw,
					  double* M)
{
  for (int n=0; n<nq*nq; n++) M[n] = 0.;
}
// [BCPenaltyJacobian]


// [penaltyData]
void Strand2dFCSPLamBc::penaltyData(const double* x,
				    double* f)
{
  double p,u,v,t,r,re;
  p    = bValue[0];
  u    = bValue[1];
  v    = bValue[2];
  t    = bValue[3];
  r    = p/(rGas*t);
  re   = p/gm1+.5*r*(u*u+v*v);
  f[0] = r;
  f[1] = r*u;
  f[2] = r*v;
  f[3] = re;
}
// [penaltyData]


// [surfaceForces]
void Strand2dFCSPLamBc::surfaceForces(const double* xs,
				      const double* ys,
				      const double* q,
				      const double* qa,
				      const double* qx,
				      const double* qy,
				      const double* qax,
				      const double* qay,
				      double* force)
{
  for (int n=0; n<ndim; n++) force[n] = 0.;
}
// [surfaceForces]


// [finalize]
void Strand2dFCSPLamBc::finalize()
{
  if (bValue) delete [] bValue;
  bValue = NULL;
}
// [finalize]
