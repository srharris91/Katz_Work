#include "Sutherland.h"

///\cond extern
extern "C" void inputtransportreadsuth_(const int&,
					const char*,
					const int&,
					double&,
					double&,
					double&,
					double&,
					double&);
///\endcond


// [Sutherland]
Sutherland::Sutherland()
{
}
// [Sutherland]


// [~Sutherland]
Sutherland::~Sutherland()
{
}
// [~Sutherland]


// [initialize]
void Sutherland::initialize(const int& comp,
			    const string& inputFile,
			    State* state0)
{
  int i=comp+1;//Fortran is 1-based
  inputtransportreadsuth_(inputFile.size(),
			  inputFile.c_str(),
			  i,
			  prn,
			  prnT,
			  t0,
			  mu0,
			  s);
  state = state0;
}
// [initialize]


// [getPrn]
void Sutherland::getPrn(double* prn0)
{
  prn0[0] = prn;
}
// [getPrn]


// [getPrnT]
void Sutherland::getPrnT(double* prnT0)
{
  prnT0[0] = prnT;
}
// [getPrnT]


// [getViscosity]
void Sutherland::getViscosity(const int& npts,
			      const double* p,
			      const double* t,
			      double* mu)
{
  double a;
  for (int n=0; n<npts; n++){
    a = sqrt(t[n]/t0);
    mu[n] = mu0*pow(a,3)*(t0+s)/(t[n]+s);
  }
}
// [getViscosity]


// [getViscosityX]
void Sutherland::getViscosityX(const int& npts,
			       const double* p,
			       const double* t,
			       const double* px,
			       const double* tx,
			       double* mux)
{
  double a = sqrt(t0);
  a = mu0*(t0+s)/pow(a,3);
  for (int n=0; n<npts; n++)
    mux[n] = a*sqrt(t[n])*tx[n]/(t[n]+s)*(1.5-t[n]/(t[n]+s));
}
// [getViscosityX]


// [getViscosityY]
void Sutherland::getViscosityY(const int& npts,
			       const double* p,
			       const double* t,
			       const double* py,
			       const double* ty,
			       double* muy)
{
  double a = sqrt(t0);
  a = mu0*(t0+s)/pow(a,3);
  for (int n=0; n<npts; n++)
    muy[n] = a*sqrt(t[n])*ty[n]/(t[n]+s)*(1.5-t[n]/(t[n]+s));
}
// [getViscosityY]


// [getViscosityZ]
void Sutherland::getViscosityZ(const int& npts,
			       const double* p,
			       const double* t,
			       const double* pz,
			       const double* tz,
			       double* muz)
{
  double a = sqrt(t0);
  a = mu0*(t0+s)/pow(a,3);
  for (int n=0; n<npts; n++)
    muz[n] = a*sqrt(t[n])*tz[n]/(t[n]+s)*(1.5-t[n]/(t[n]+s));
}
// [getViscosityZ]


// [getDViscosityDT]
void Sutherland::getDViscosityDT(const int& npts,
				 const double* p,
				 const double* t,
				 double* dmudT)
{
  double a;
  for (int n=0; n<npts; n++){
    a        = sqrt(t[n]/t0);
    dmudT[n] = mu0*(1.5/t0*a*(t0+s)/(t[n]+s)
		    -pow(a,3)*(t0+s)/((t[n]+s)*(t[n]+s)));
  }
}
// [getDViscosityDT]


// [getConductivity]
void Sutherland::getConductivity(const int& npts,
				 const double* p,
				 const double* t,
				 double* k)
{
  double* Cp = new double[npts];
  state->getCp(npts,
	       p,
	       t,
	       Cp);
  double b,a = mu0*(t0+s)/prn;
  for (int n=0; n<npts; n++){
    b = sqrt(t[n]/t0);
    k[n] = a*Cp[n]*pow(b,3)/(t[n]+s);
  }
  delete [] Cp;
}
// [getConductivity]


// [getConductivityX]
void Sutherland::getConductivityX(const int& npts,
				  const double* p,
				  const double* t,
				  const double* px,
				  const double* tx,
				  double* kx)
{
  double* Cp = new double[npts];
  state->getCp(npts,
	       p,
	       t,
	       Cp);
  double b = sqrt(t0);
  double a = mu0*(t0+s)/(prn*pow(b,3));
  for (int n=0; n<npts; n++)
    kx[n] = a*Cp[n]*sqrt(t[n])*tx[n]/(t[n]+s)*(1.5-t[n]/(t[n]+s));

  delete [] Cp;
}
// [getConductivityX]


// [getConductivityY]
void Sutherland::getConductivityY(const int& npts,
				  const double* p,
				  const double* t,
				  const double* py,
				  const double* ty,
				  double* ky)
{
  double* Cp = new double[npts];
  state->getCp(npts,
	       p,
	       t,
	       Cp);
  double b = sqrt(t0);
  double a = mu0*(t0+s)/(prn*pow(b,3));
  for (int n=0; n<npts; n++)
    ky[n] = a*Cp[n]*sqrt(t[n])*ty[n]/(t[n]+s)*(1.5-t[n]/(t[n]+s));

  delete [] Cp;
}
// [getConductivityY]


// [getConductivityZ]
void Sutherland::getConductivityZ(const int& npts,
				  const double* p,
				  const double* t,
				  const double* pz,
				  const double* tz,
				  double* kz)
{
  double* Cp = new double[npts];
  state->getCp(npts,
	       p,
	       t,
	       Cp);
  double b = sqrt(t0);
  double a = mu0*(t0+s)/(prn*pow(b,3));
  for (int n=0; n<npts; n++)
    kz[n] = a*Cp[n]*sqrt(t[n])*tz[n]/(t[n]+s)*(1.5-t[n]/(t[n]+s));

  delete [] Cp;
}
// [getConductivityZ]


// [getDConductivityDT]
void Sutherland::getDConductivityDT(const int& npts,
				    const double* p,
				    const double* t,
				    double* dkdT)
{
  double* Cp = new double[npts];
  state->getCp(npts,
	       p,
	       t,
	       Cp);
  double a;
  for (int n=0; n<npts; n++){
    a       = sqrt(t[n]/t0);
    dkdT[n] = Cp[n]*mu0/prn*(1.5/t0*a*(t0+s)/(t[n]+s)
			     -pow(a,3)*(t0+s)/((t[n]+s)*(t[n]+s)));
  }
  delete [] Cp;
}
// [getDConductivityDT]


// [finalize]
void Sutherland::finalize()
{
}
// [finalize]
