#include "Incompressible.h"

///\cond extern
extern "C"
{
  void inputstatereadincompressible_(const int&,
				     const char*,
				     const int&,
				     double&,
				     double&,
				     double&);
}
///\endcond


// [Incompressible]
Incompressible::Incompressible()
{
  rho = 0.;
  Cp  = 0.;
  Cv  = 0.;
  c   = 1.e6; //something very large
}
// [Incompressible]


// [~Incompressible]
Incompressible::~Incompressible()
{
}
// [~Incompressible]


// [initialize]
void Incompressible::initialize(const int& comp,
				const string& inputFile)
{
  int i=comp+1;//Fortran is 1-based
  inputstatereadincompressible_(inputFile.size(),
				inputFile.c_str(),
				i,
				rho,
				Cp,
				Cv);
}
// [initialize]


// [getGamma]
void Incompressible::getGamma(const int& npts,
			      const double* p,
			      const double* t,
			      double* gamma0)
{
  for (int n=0; n<npts; n++) gamma0[n] = Cp/Cv;
}
// [getGamma]


// [getCp]
void Incompressible::getCp(const int& npts,
			   const double* p,
			   const double* t,
			   double* Cpn)
{
  for (int n=0; n<npts; n++) Cpn[n] = Cp;
}
// [getCp]


// [getSoundSpeed]
void Incompressible::getSoundSpeed(const int& npts,
				   const double* p,
				   const double* t,
				   double* cn)
{
  for (int n=0; n<npts; n++) cn[n] = c;
}
// [getSoundSpeed]


// [getRhoP]
void Incompressible::getRhoP(const int& npts,
			     const double* p,
			     const double* t,
			     double* rp)
{
  for (int n=0; n<npts; n++) rp[n] = 0.;
}
// [getRhoP]


// [getRhoT]
void Incompressible::getRhoT(const int& npts,
			     const double* p,
			     const double* t,
			     double* rt)
{
  for (int n=0; n<npts; n++) rt[n] = 0.;
}
// [getRhoT]


// [getHP]
void Incompressible::getHP(const int& npts,
			   const double* p,
			   const double* t,
			   double* hp)
{
  for (int n=0; n<npts; n++) hp[n] = 0.;
}
// [getHP]


// [getHT]
void Incompressible::getHT(const int& npts,
			   const double* p,
			   const double* t,
			   double* ht)
{
  for (int n=0; n<npts; n++) ht[n] = Cp;
}
// [getHT]


// [getDensity]
void Incompressible::getDensity(const int& npts,
				const double* p,
				const double* t,
				double* rhon)
{
  for (int n=0; n<npts; n++) rhon[n] = rho;
}
// [getDensity]


// [getDensityX]
void Incompressible::getDensityX(const int& npts,
				 const double* p,
				 const double* t,
				 const double* px,
				 const double* tx,
				 double* rhox)
{
  for (int n=0; n<npts; n++) rhox[n] = 0.;
}
// [getDensityX]


// [getDensityY]
void Incompressible::getDensityY(const int& npts,
				 const double* p,
				 const double* t,
				 const double* py,
				 const double* ty,
				 double* rhoy)
{
  for (int n=0; n<npts; n++) rhoy[n] = 0.;
}
// [getDensityY]


// [getDensityZ]
void Incompressible::getDensityZ(const int& npts,
				 const double* p,
				 const double* t,
				 const double* pz,
				 const double* tz,
				 double* rhoz)
{
  for (int n=0; n<npts; n++) rhoz[n] = 0.;
}
// [getDensityZ]


// [getEnthalpy]
void Incompressible::getEnthalpy(const int& npts,
				 const double* p,
				 const double* t,
				 double* h)
{
  for (int n=0; n<npts; n++) h[n] = Cp*t[n];
}
// [getEnthalpy]


// [getEnthalpyX]
void Incompressible::getEnthalpyX(const int& npts,
				  const double* p,
				  const double* t,
				  const double* px,
				  const double* tx,
				  double* hx)
{
  for (int n=0; n<npts; n++) hx[n] = Cp*tx[n];
}
// [getEnthalpyX]


// [getEnthalpyY]
void Incompressible::getEnthalpyY(const int& npts,
				  const double* p,
				  const double* t,
				  const double* py,
				  const double* ty,
				  double* hy)
{
  for (int n=0; n<npts; n++) hy[n] = Cp*ty[n];
}
// [getEnthalpyY]


// [getEnthalpyZ]
void Incompressible::getEnthalpyZ(const int& npts,
				  const double* p,
				  const double* t,
				  const double* pz,
				  const double* tz,
				  double* hz)
{
  for (int n=0; n<npts; n++) hz[n] = Cp*tz[n];
}
// [getEnthalpyZ]


// [finalize]
void Incompressible::finalize()
{
}
// [finalize]
