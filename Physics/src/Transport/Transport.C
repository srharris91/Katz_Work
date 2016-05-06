#include "Transport.h"


// [Transport]
Transport::Transport()
{
}
// [Transport]


// [~Transport]
Transport::~Transport()
{
}
// [~Transport]


// [initialize]
void Transport::initialize(const int& ncomp0,
			   const int* itransport,
			   const string& inputFile,
			   State* state)
{
  ncomp = ncomp0;
  transportType = new TransportType*[ncomp];
  for (int n=0; n<ncomp; n++){
    if (itransport[n] == 1) transportType[n] = new Sutherland;
    else{
      cout << "\nTransport type not recognized in initialize. Terminating."
	   << endl;
      exit(0);
    }
    transportType[n]->initialize(n,
				 inputFile,
				 state);
  }
}
// [initialize]


// [getPrn]
void Transport::getPrn(double* Prn)
{
  for (int n=0; n<ncomp; n++) transportType[n]->getPrn(Prn+n);
}
// [getPrn]


// [getPrnT]
void Transport::getPrnT(double* PrnT)
{
  for (int n=0; n<ncomp; n++) transportType[n]->getPrnT(PrnT+n);
}
// [getPrnT]


// [getViscosity]
void Transport::getViscosity(const int& npts,
			     const double* p,
			     const double* t,
			     double* mu)
{
  for (int n=0; n<ncomp; n++)
    transportType[n]->getViscosity(npts,
				   p +n*npts,
				   t +n*npts,
				   mu+n*npts);
}
// [getViscosity]


// [getViscosityX]
void Transport::getViscosityX(const int& npts,
			      const double* p,
			      const double* t,
			      const double* px,
			      const double* tx,
			      double* mux)
{
  for (int n=0; n<ncomp; n++)
    transportType[n]->getViscosityX(npts,
				    p  +n*npts,
				    t  +n*npts,
				    px +n*npts,
				    tx +n*npts,
				    mux+n*npts);
}
// [getViscosityX]


// [getViscosityY]
void Transport::getViscosityY(const int& npts,
			      const double* p,
			      const double* t,
			      const double* py,
			      const double* ty,
			      double* muy)
{
  for (int n=0; n<ncomp; n++)
    transportType[n]->getViscosityY(npts,
				    p  +n*npts,
				    t  +n*npts,
				    py +n*npts,
				    ty +n*npts,
				    muy+n*npts);
}
// [getViscosityY]


// [getViscosityZ]
void Transport::getViscosityZ(const int& npts,
			      const double* p,
			      const double* t,
			      const double* pz,
			      const double* tz,
			      double* muz)
{
  for (int n=0; n<ncomp; n++)
    transportType[n]->getViscosityZ(npts,
				    p  +n*npts,
				    t  +n*npts,
				    pz +n*npts,
				    tz +n*npts,
				    muz+n*npts);
}
// [getViscosityZ]


// [getDViscosityDT]
void Transport::getDViscosityDT(const int& npts,
				const double* p,
				const double* t,
				double* dmudT)
{
  for (int n=0; n<ncomp; n++)
    transportType[n]->getDViscosityDT(npts,
				      p    +n*npts,
				      t    +n*npts,
				      dmudT+n*npts);
}
// [getDViscosityDT]


// [getConductivity]
void Transport::getConductivity(const int& npts,
				const double* p,
				const double* t,
				double* k)
{
  for (int n=0; n<ncomp; n++)
    transportType[n]->getConductivity(npts,
				      p +n*npts,
				      t +n*npts,
				      k+n*npts);
}
// [getConductivity]


// [getConductivityX]
void Transport::getConductivityX(const int& npts,
				 const double* p,
				 const double* t,
				 const double* px,
				 const double* tx,
				 double* kx)
{
  for (int n=0; n<ncomp; n++)
    transportType[n]->getConductivityX(npts,
				       p  +n*npts,
				       t  +n*npts,
				       px +n*npts,
				       tx +n*npts,
				       kx+n*npts);
}
// [getConductivityX]


// [getConductivityY]
void Transport::getConductivityY(const int& npts,
				 const double* p,
				 const double* t,
				 const double* py,
				 const double* ty,
				 double* ky)
{
  for (int n=0; n<ncomp; n++)
    transportType[n]->getConductivityY(npts,
				       p  +n*npts,
				       t  +n*npts,
				       py +n*npts,
				       ty +n*npts,
				       ky+n*npts);
}
// [getConductivityY]


// [getConductivityZ]
void Transport::getConductivityZ(const int& npts,
				 const double* p,
				 const double* t,
				 const double* pz,
				 const double* tz,
				 double* kz)
{
  for (int n=0; n<ncomp; n++)
    transportType[n]->getConductivityZ(npts,
				       p  +n*npts,
				       t  +n*npts,
				       pz +n*npts,
				       tz +n*npts,
				       kz+n*npts);
}
// [getConductivityZ]


// [getDConductivityDT]
void Transport::getDConductivityDT(const int& npts,
				   const double* p,
				   const double* t,
				   double* dkdT)
{
  for (int n=0; n<ncomp; n++)
    transportType[n]->getDConductivityDT(npts,
					 p   +n*npts,
					 t   +n*npts,
					 dkdT+n*npts);
}
// [getDConductivityDT]


// [finalize]
void Transport::finalize()
{
  if (transportType){
    for (int n=0; n<ncomp; n++){
      transportType[n]->finalize();
      delete transportType[n];
      transportType[n] = NULL;
    }
    delete [] transportType;
    transportType = NULL;
  }
}
// [finalize]
