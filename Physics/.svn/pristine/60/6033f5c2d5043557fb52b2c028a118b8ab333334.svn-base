#include "State.h"


// [State]
State::State()
{
}
// [State]


// [~State]
State::~State()
{
}
// [~State]


// [initialize]
void State::initialize(const int& ncomp0,
		       const int* istate,
		       const string& inputFile)
{
  ncomp = ncomp0;
  stateType = new StateType*[ncomp];
  for (int n=0; n<ncomp; n++){
    if      (istate[n] == 1) stateType[n] = new IdealGas;
    else if (istate[n] == 2) stateType[n] = new Incompressible;
    else{
      cout << "\nState type not recognized in initialize. Terminating."
	   << endl;
      exit(0);
    }
    stateType[n]->initialize(n,inputFile);
  }
}
// [initialize]


// [getGamma]
void State::getGamma(const int& npts,
		     const double* p,
		     const double* t,
		     double* gamma)
{
  for (int n=0; n<ncomp; n++)
    stateType[n]->getGamma(npts,
			   p    +n*npts,
			   t    +n*npts,
			   gamma+n*npts);
}
// [getGamma]


// [getGamma2]
void State::getGamma(double* gamma)
{
  for (int n=0; n<ncomp; n++) stateType[n]->getGamma(gamma+n);
}
// [getGamma2]


// [getRGas]
void State::getRGas(double* rGas)
{
  for (int n=0; n<ncomp; n++) stateType[n]->getRGas(rGas+n);
}
// [getRGas]


// [getCp]
void State::getCp(const int& npts,
		  const double* p,
		  const double* t,
		  double* Cp)
{
  for (int n=0; n<ncomp; n++)
    stateType[n]->getCp(npts,
			p +n*npts,
			t +n*npts,
			Cp+n*npts);
}
// [getCp]


// [getSoundSpeed]
void State::getSoundSpeed(const int& npts,
			  const double* p,
			  const double* t,
			  double* c)
{
  for (int n=0; n<ncomp; n++)
    stateType[n]->getSoundSpeed(npts,
				p+n*npts,
				t+n*npts,
				c+n*npts);
}
// [getSoundSpeed]


// [getRhoP]
void State::getRhoP(const int& npts,
		    const double* p,
		    const double* t,
		    double* rp)
{
  for (int n=0; n<ncomp; n++)
    stateType[n]->getRhoP(npts,
			  p +n*npts,
			  t +n*npts,
			  rp+n*npts);
}
// [getRhoP]


// [getRhoT]
void State::getRhoT(const int& npts,
		    const double* p,
		    const double* t,
		    double* rt)
{
  for (int n=0; n<ncomp; n++)
    stateType[n]->getRhoT(npts,
			  p +n*npts,
			  t +n*npts,
			  rt+n*npts);
}
// [getRhoT]


// [getHP]
void State::getHP(const int& npts,
		  const double* p,
		  const double* t,
		  double* hp)
{
  for (int n=0; n<ncomp; n++)
    stateType[n]->getHP(npts,
			p +n*npts,
			t +n*npts,
			hp+n*npts);
}
// [getHP]


// [getHT]
void State::getHT(const int& npts,
		  const double* p,
		  const double* t,
		  double* ht)
{
  for (int n=0; n<ncomp; n++)
    stateType[n]->getHT(npts,
			p +n*npts,
			t +n*npts,
			ht+n*npts);
}
// [getHT]


// [getDensity]
void State::getDensity(const int& npts,
		       const double* p,
		       const double* t,
		       double* rho)
{
  for (int n=0; n<ncomp; n++)
    stateType[n]->getDensity(npts,
			     p  +n*npts,
			     t  +n*npts,
			     rho+n*npts);
}
// [getDensity]


// [getDensityX]
void State::getDensityX(const int& npts,
			const double* p,
			const double* t,
			const double* px,
			const double* tx,
			double* rhox)
{
  for (int n=0; n<ncomp; n++)
    stateType[n]->getDensityX(npts,
			      p   +n*npts,
			      t   +n*npts,
			      px  +n*npts,
			      tx  +n*npts,
			      rhox+n*npts);
}
// [getDensityX]


// [getDensityY]
void State::getDensityY(const int& npts,
			const double* p,
			const double* t,
			const double* py,
			const double* ty,
			double* rhoy)
{
  for (int n=0; n<ncomp; n++)
    stateType[n]->getDensityY(npts,
			      p   +n*npts,
			      t   +n*npts,
			      py  +n*npts,
			      ty  +n*npts,
			      rhoy+n*npts);
}
// [getDensityY]


// [getDensityZ]
void State::getDensityZ(const int& npts,
			const double* p,
			const double* t,
			const double* pz,
			const double* tz,
			double* rhoz)
{
  for (int n=0; n<ncomp; n++)
    stateType[n]->getDensityZ(npts,
			      p   +n*npts,
			      t   +n*npts,
			      pz  +n*npts,
			      tz  +n*npts,
			      rhoz+n*npts);
}
// [getDensityZ]


// [getEnthalpy]
void State::getEnthalpy(const int& npts,
			const double* p,
			const double* t,
			double* h)
{
  for (int n=0; n<ncomp; n++)
    stateType[n]->getEnthalpy(npts,
			      p+n*npts,
			      t+n*npts,
			      h+n*npts);
}
// [getEnthalpy]


// [getEnthalpyMix]
void State::getEnthalpy(const int& npts,
			const double* g,
			const double* p,
			const double* t,
			double* h)
{
  double* rhoK = new double[npts*ncomp];
  double* hK   = new double[npts*ncomp];
  for (int n=0; n<ncomp; n++){
    stateType[n]->getEnthalpy(npts,
			      p,
			      t,
			      hK+n*npts);
  }
  
  for (int n=0; n<npts; n++) h[n] = 0.;
  
  int i;
  for (int n=0; n<npts; n++){
  for (int m=0; m<ncomp; m++){
    i     = m*npts+n;
    h[n] += g[i]*hK[i];
  }
  }
}
// [getEnthalpyMix]


// [getEnthalpyX]
void State::getEnthalpyX(const int& npts,
			 const double* p,
			 const double* t,
			 const double* px,
			 const double* tx,
			 double* hx)
{
  for (int n=0; n<ncomp; n++)
    stateType[n]->getEnthalpyX(npts,
			      p +n*npts,
			      t +n*npts,
			      px+n*npts,
			      tx+n*npts,
			      hx+n*npts);
}
// [getEnthalpyX]


// [getEnthalpyY]
void State::getEnthalpyY(const int& npts,
			 const double* p,
			 const double* t,
			 const double* py,
			 const double* ty,
			 double* hy)
{
  for (int n=0; n<ncomp; n++)
    stateType[n]->getEnthalpyY(npts,
			      p +n*npts,
			      t +n*npts,
			      py+n*npts,
			      ty+n*npts,
			      hy+n*npts);
}
// [getEnthalpyY]


// [getEnthalpyZ]
void State::getEnthalpyZ(const int& npts,
			 const double* p,
			 const double* t,
			 const double* pz,
			 const double* tz,
			 double* hz)
{
  for (int n=0; n<ncomp; n++)
    stateType[n]->getEnthalpyZ(npts,
			      p +n*npts,
			      t +n*npts,
			      pz+n*npts,
			      tz+n*npts,
			      hz+n*npts);
}
// [getEnthalpyZ]


// [getDensityMix]
void State::getDensity(const int& npts,
			     const double* g,
			     const double* p,
			     const double* t,
			     double* rho)
{
  double* rhoK = new double[npts*ncomp];
  for (int n=0; n<ncomp; n++){
    stateType[n]->getDensity(npts,
			     p,
			     t,
			     rhoK+n*npts);
  }

  for (int n=0; n<npts; n++) rho[n] = 0.;

  int i;
  for (int n=0; n<npts; n++){
  for (int m=0; m<ncomp; m++){
    i       = m*npts+n;
    rho[n] += g[i]/rhoK[i];
  }
  }
  for (int n=0; n<npts; n++)rho[n] = 1./rho[n];
}
// [getDensityMix]


// [finalize]
void State::finalize()
{
  if (stateType){
    for (int n=0; n<ncomp; n++){
      stateType[n]->finalize();
      delete stateType[n];
      stateType[n] = NULL;
    }
    delete [] stateType;
    stateType = NULL;
  }
}
// [finalize]
