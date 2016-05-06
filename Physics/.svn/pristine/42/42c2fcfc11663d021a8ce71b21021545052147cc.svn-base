#include "Mms2d.h"

///\cond extern
extern "C"
{
  void inputsolutionreadmms2d_(const int&,
			       const char*,
			       const int&,
			       double&,
			       double&,
			       double*);
}
///\endcond


// [Mms2d]
Mms2d::Mms2d()
{
}
// [Mms2d]


// [~Mms2d]
Mms2d::~Mms2d()
{
}
// [~Mms2d]


// [initialize]
void Mms2d::initialize(const int& nq0,
		       const int& ndim0,
		       const string& inputFile,
		       State* state0,
		       Transport* transport0)
{
  nq = nq0;
  ndim = ndim0;
  rValue = new double[nq];
  inputsolutionreadmms2d_(inputFile.size(),
			  inputFile.c_str(),
			  nq,
			  period,
			  amplitude,
			  rValue);
  state     = state0;
  transport = transport0;

  ax  = new double[nq];
  ay  = new double[nq];
  axy = new double[nq];
  bx  = new double[nq];
  by  = new double[nq];
  bxy = new double[nq];
  cx  = new double[nq];
  cy  = new double[nq];
  cxy = new double[nq];

  double pi = 4.*atan(1.);
  double A = amplitude;
  double a = amplitude/3.;
  double B = 2.*pi/period;
  double b,c;
  srand(1);
  // generate random number MMS coefficients (since the seed is the same
  // for all instances, the same coefficients will be generated each time).
  for (int n=0; n<nq; n++){

    b      = double(rand())/double(RAND_MAX);
    //b = .05/double(n+1);
    b      = 2.*a*(b-.5); // random number between -a and a
    ax[n]  = b;
    b      = double(rand())/double(RAND_MAX);
    //b = .15/double(n+1);
    b      = 2.*a*(b-.5); // random number between -a and a
    ay[n]  = b;
    b      = double(rand())/double(RAND_MAX);
    //b = .25/double(n+1);
    b      = 2.*a*(b-.5); // random number between -a and a
    axy[n] = b;

    //ax[n]  = 0.;
    //ay[n]  = 0.;
    //axy[n] = 0.;

    b      = double(rand())/double(RAND_MAX);
    //b = .35/double(n+1);
    c      = 0.;
    if (b < .5) c = 1.; // 0. or 1. randomly
    bx[n]  = c*.5*pi;
    b      = double(rand())/double(RAND_MAX);
    //b = .45/double(n+1);
    c      = 0.;
    if (b < .5) c = 1.; // 0. or 1. randomly
    by[n]  = c*.5*pi;
    b      = double(rand())/double(RAND_MAX);
    //b = .55/double(n+1);
    c      = 0.;
    if (b < .5) c = 1.; // 0. or 1. randomly
    bxy[n] = c*.5*pi;

    b      = double(rand())/double(RAND_MAX);
    //b = .65/double(n+1);
    b      = 1.+2.*A*(b-.5); // random number between 1.-A and 1.+A
    cx[n]  = b*B;
    b      = double(rand())/double(RAND_MAX);
    //b = .75/double(n+1);
    b      = 1.+2.*A*(b-.5); // random number between 1.-A and 1.+A
    cy[n]  = b*B;
    b      = double(rand())/double(RAND_MAX);
    //b = .85/double(n+1);
    b      = 1.+2.*A*(b-.5); // random number between 1.-A and 1.+A
    cxy[n] = b*B/period;
  }
}
void Mms2d::initialize(const int& nq0,
		       const int& ndim0,
		       const string& inputFile)
{
  nq = nq0;
  ndim = ndim0;
  rValue = new double[nq];
  inputsolutionreadmms2d_(inputFile.size(),
			  inputFile.c_str(),
			  nq,
			  period,
			  amplitude,
			  rValue);

  ax  = new double[nq];
  ay  = new double[nq];
  axy = new double[nq];
  bx  = new double[nq];
  by  = new double[nq];
  bxy = new double[nq];
  cx  = new double[nq];
  cy  = new double[nq];
  cxy = new double[nq];

  double pi = 4.*atan(1.);
  double A = amplitude;
  double a = amplitude/3.;
  double B = 2.*pi/period;
  double b,c;
  srand(1);
  // generate random number MMS coefficients (since the seed is the same
  // for all instances, the same coefficients will be generated each time).
  for (int n=0; n<nq; n++){

    b      = double(rand())/double(RAND_MAX);
    //b = .05/double(n+1);
    b      = 2.*a*(b-.5); // random number between -a and a
    ax[n]  = b;
    b      = double(rand())/double(RAND_MAX);
    //b = .15/double(n+1);
    b      = 2.*a*(b-.5); // random number between -a and a
    ay[n]  = b;
    b      = double(rand())/double(RAND_MAX);
    //b = .25/double(n+1);
    b      = 2.*a*(b-.5); // random number between -a and a
    axy[n] = b;

    //ax[n]  = 0.;
    //ay[n]  = 0.;
    //axy[n] = 0.;

    b      = double(rand())/double(RAND_MAX);
    //b = .35/double(n+1);
    c      = 0.;
    if (b < .5) c = 1.; // 0. or 1. randomly
    bx[n]  = c*.5*pi;
    b      = double(rand())/double(RAND_MAX);
    //b = .45/double(n+1);
    c      = 0.;
    if (b < .5) c = 1.; // 0. or 1. randomly
    by[n]  = c*.5*pi;
    b      = double(rand())/double(RAND_MAX);
    //b = .55/double(n+1);
    c      = 0.;
    if (b < .5) c = 1.; // 0. or 1. randomly
    bxy[n] = c*.5*pi;

    b      = double(rand())/double(RAND_MAX);
    //b = .65/double(n+1);
    b      = 1.+2.*A*(b-.5); // random number between 1.-A and 1.+A
    cx[n]  = b*B;
    b      = double(rand())/double(RAND_MAX);
    //b = .75/double(n+1);
    b      = 1.+2.*A*(b-.5); // random number between 1.-A and 1.+A
    cy[n]  = b*B;
    b      = double(rand())/double(RAND_MAX);
    //b = .85/double(n+1);
    b      = 1.+2.*A*(b-.5); // random number between 1.-A and 1.+A
    cxy[n] = b*B/period;
  }
}
// [initialize]


// [getQ]
void Mms2d::getQ(const int& npts,
		 const double* xy,
		 double* q)
{
  int i,iq,ix;
  double x,y;
  for (int n=0; n<npts; n++){
    iq     = n*nq;
    ix     = n*ndim;
    x      = xy[ix  ];
    y      = xy[ix+1];
    for (int k=0; k<nq; k++){
      q[iq+k] = rValue[k]*(1.+ax [k]*sin(bx [k]-cx [k]*x)
			     +ay [k]*sin(by [k]-cy [k]*y)
			     +axy[k]*sin(bxy[k]-cxy[k]*x*y));
    }}
}
// [getQ]


// [getQx]
void Mms2d::getQx(const int& npts,
		  const double* xy,
		  double* qx)
{
  int i,iq,ix;
  double x,y;
  for (int n=0; n<npts; n++){
    iq     = n*nq;
    ix     = n*ndim;
    x      = xy[ix  ];
    y      = xy[ix+1];
    for (int k=0; k<nq; k++)
      qx[iq+k] =-rValue[k]*(ax [k]*cx [k]  *cos(bx [k]-cx [k]*x)
			   +axy[k]*cxy[k]*y*cos(bxy[k]-cxy[k]*x*y));
  }
}
// [getQx]


// [getQy]
void Mms2d::getQy(const int& npts,
		  const double* xy,
		  double* qy)
{
  int i,iq,ix;
  double x,y;
  for (int n=0; n<npts; n++){
    iq     = n*nq;
    ix     = n*ndim;
    x      = xy[ix  ];
    y      = xy[ix+1];
    for (int k=0; k<nq; k++)
      qy[iq+k] =-rValue[k]*(ay [k]*cy [k]  *cos(by [k]-cy [k]*y)
			   +axy[k]*cxy[k]*x*cos(bxy[k]-cxy[k]*x*y));
  }
}
// [getQy]


// [getQxx]
void Mms2d::getQxx(const int& npts,
		   const double* xy,
		   double* qxx)
{
  int i,iq,ix;
  double x,y;
  for (int n=0; n<npts; n++){
    iq     = n*nq;
    ix     = n*ndim;
    x      = xy[ix  ];
    y      = xy[ix+1];
    for (int k=0; k<nq; k++)
      qxx[iq+k] =-rValue[k]*(ax [k]*cx [k]*cx [k]    *sin(bx [k]-cx [k]*x)
			    +axy[k]*cxy[k]*cxy[k]*y*y*sin(bxy[k]-cxy[k]*x*y));
  }
}
// [getQxx]


// [getQyy]
void Mms2d::getQyy(const int& npts,
		   const double* xy,
		   double* qyy)
{
  int i,iq,ix;
  double x,y;
  for (int n=0; n<npts; n++){
    iq     = n*nq;
    ix     = n*ndim;
    x      = xy[ix  ];
    y      = xy[ix+1];
    for (int k=0; k<nq; k++)
      qyy[iq+k] =-rValue[k]*(ay [k]*cy [k]*cy [k]    *sin(by [k]-cy [k]*y)
			    +axy[k]*cxy[k]*cxy[k]*x*x*sin(bxy[k]-cxy[k]*x*y));
  }
}
// [getQyy]


// [getQxy]
void Mms2d::getQxy(const int& npts,
		   const double* xy,
		   double* qxy)
{
  int i,iq,ix;
  double x,y;
  for (int n=0; n<npts; n++){
    iq     = n*nq;
    ix     = n*ndim;
    x      = xy[ix  ];
    y      = xy[ix+1];
    for (int k=0; k<nq; k++)
      qxy[iq+k] =-rValue[k]*(axy[k]*cxy[k]*(cxy[k]*x*y*sin(bxy[k]-cxy[k]*x*y)+
					    cos(bxy[k]-cxy[k]*x*y)));
  }
}
// [getQxy]


// [finalize]
void Mms2d::finalize()
{
  if (ax ) delete [] ax;
  if (ay ) delete [] ay;
  if (axy) delete [] axy;
  if (bx ) delete [] bx;
  if (by ) delete [] by;
  if (bxy) delete [] bxy;
  if (cx ) delete [] cx;
  if (cy ) delete [] cy;
  if (cxy) delete [] cxy;
  ax  = NULL;
  ay  = NULL;
  axy = NULL;
  bx  = NULL;
  by  = NULL;
  bxy = NULL;
  cx  = NULL;
  cy  = NULL;
  cxy = NULL;
}
// [finalize]
