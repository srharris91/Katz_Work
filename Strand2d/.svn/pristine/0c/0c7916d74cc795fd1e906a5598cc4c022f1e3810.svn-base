#include "StrandSPLam.h"


void StrandSPLam::rhsDisFluxCoarse(const int& npts,
				   const double* A,
				   const double* xv,
				   const double* ql,
				   const double* qr,
				   double* f)
{
  // scalar dissipative flux terms
  int iq,iA;
  double Ax,Ay,rhl,rhr,rl,rr,dd,ua,va,ha,cs,qs,pl,pr,rul,rvl,re,
    rqq,rur,rvr,cc,qq,sr,dlim=.25;
  for (int n=0; n<npts; n++){
    iq      = nq  *n;
    iA      = ndim*n;
    Ax      = A[iA  ];
    Ay      = A[iA+1];

    rl      = ql[iq  ];
    rul     = ql[iq+1];
    rvl     = ql[iq+2];
    re      = ql[iq+3];
    rqq     =(rul*rul+rvl*rvl)/rl;
    pl      = gm1*(re-.5*rqq);
    rhl     = re+pl;
    rr      = qr[iq  ];
    rur     = qr[iq+1];
    rvr     = qr[iq+2];
    re      = qr[iq+3];
    rqq     =(rur*rur+rvr*rvr)/rr;
    pr      = gm1*(re-.5*rqq);
    rhr     = re+pr;

    rl      = sqrt(rl);
    rr      = sqrt(rr);
    dd      = 1./(rl+rr);
    rl      = 1./rl;
    rr      = 1./rr;
    ua      =(rul*rl+rur*rr)*dd;
    va      =(rvl*rl+rvr*rr)*dd;
    ha      =(rhl*rl+rhr*rr)*dd;
    qq      = .5*(ua*ua+va*va);
    cc      = gm1*(ha-qq);

    cs      = sqrt(cc*(Ax*Ax+Ay*Ay));
    qs      = Ax*ua+Ay*va-xv[n];
    sr      = dlim*(fabs(qs)+cs);
 
    f[iq  ] = sr*(qr[iq  ]-ql[iq  ]);
    f[iq+1] = sr*(qr[iq+1]-ql[iq+1]);
    f[iq+2] = sr*(qr[iq+2]-ql[iq+2]);
    f[iq+3] = sr*(qr[iq+3]-ql[iq+3]);
  }
}
