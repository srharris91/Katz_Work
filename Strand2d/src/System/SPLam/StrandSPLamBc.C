#include "StrandSPLamBc.h"


// [StrandSPLamBc]
StrandSPLamBc::StrandSPLamBc()
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
  transport = NULL;
  solution = NULL;
}
// [StrandSPLamBc]


// [~StrandSPLamBc]
StrandSPLamBc::~StrandSPLamBc()
{
}
// [~StrandSPLamBc]


// [initialize]
void StrandSPLamBc::initialize(const int& nq0,
			       const int& nqa0,
			       const int& ndim0,
			       const int& ncomp0,
			       const double& gamma0,
			       const double& rGas0,
			       double* bValue0,
			       Transport* transport0,
			       Solution* solution0)
{
  nq = nq0;
  nqa = nqa0;
  ndim = ndim0;
  ncomp = ncomp0;
  gamma = gamma0;
  rGas = rGas0;
  gm1 = gamma-1.;
  ggm1 = gamma/gm1;
  if (!bValue) bValue = new double[nq];
  for (int n=0; n<nq; n++) bValue[n] = bValue0[n];
  transport = transport0;
  solution = solution0;
}
// [initialize]


// [BCVector]
void StrandSPLamBc::BCVector(const double* nx,
			     const double* wx,
			     const double* qe,
			     const double* qae,
			     const double* q,
			     const double* qa,
			     double* r)
{
  for (int n=0; n<nq; n++) r[n] = 0;
}
// [BCVector]


// [BCVectorSelfJacobian]
void StrandSPLamBc::BCVectorSelfJacobian(const double* nx,
					 const double* qe,
					 const double* qae,
					 const double* q,
					 const double* qa,
					 double* M)
{
  M[0 ] = 1.;
  M[1 ] = 0.;
  M[2 ] = 0.;
  M[3 ] = 0.;
  M[4 ] = 0.;
  M[5 ] = 1.;
  M[6 ] = 0.;
  M[7 ] = 0.;
  M[8 ] = 0.;
  M[9 ] = 0.;
  M[10] = 1.;
  M[11] = 0.;
  M[12] = 0.;
  M[13] = 0.;
  M[14] = 0.;
  M[15] = 1.;
}
// [BCVectorSelfJacobian]


// [BCVectorInteriorJacobian]
void StrandSPLamBc::BCVectorInteriorJacobian(const double* nx,
					     const double* qe,
					     const double* qae,
					     const double* q,
					     const double* qa,
					     double* M)
{
  for (int n=0; n<nq*nq; n++) M[n] = 0.;
}
// [BCVectorInteriorJacobian]


// [surfaceForces]
void StrandSPLamBc::surfaceForces(const double* A,
				  const double* q,
				  const double* qa,
				  const double* qx,
				  const double* qax,
				  double* force)
{
  for (int n=0; n<ndim; n++) force[n] = 0.;
}
// [surfaceForces]


// [finalize]
void StrandSPLamBc::finalize()
{
  if (bValue) delete [] bValue;
  bValue = NULL;
}
// [finalize]
