#include "Strand2dFCSPTurbSA.h"


void Strand2dFCSPTurbSA::lhsVisFluxNCoeffJacobian(const int& npts,
					       const double* jac,
					       const double* xs,
					       const double* ys,
					       const double* xn,
					       const double* yn,
					       const double* q,
					       const double* qa,
					       const double* qi,
					       const double* qai,
					       double* A)
{
  int iq,iqa,iA;
  double rho,r,e,rnu,p,u,v,t,nu,mu,ku,ui,vi,Ti,nui,a,b,c,d,f,dmudT,dkdT,
    dTdq[nq],dmudq[nq],chi,chi3,fv1,mut,kt,fn,dfv1dchi,dchidq[nq],dfv1dq[nq],
    dmutdq[nq],dfndchi,dfndq[nq],dwdq[nq],dktdq[nq],dumudq[nq],dvmudq[nq],
    dkdq[nq];

  for (int n=0; n<npts; n++){
    iq       = n*nq;
    iqa      = n*nqa;
    iA       = n*nq*nq;
    rho      = q  [iq   ];
    r        = 1./rho;
    e        = q  [iq +3]*r;
    rnu      = q  [iq +4];
    p        = qa [iqa  ];
    u        = qa [iqa+1];
    v        = qa [iqa+2];
    t        = qa [iqa+3];
    nu       = qa [iqa+4];
    mu       = qa [iqa+5];
    ku       = qa [iqa+6];
    ui       = qai[iqa+1];
    vi       = qai[iqa+2];
    Ti       = qai[iqa+3];
    nui      = qai[iqa+4];

    a        = xs[n]*xs[n]+ys[n]*ys[n];
    b        =(a+ys[n]*ys[n]/3.)/jac[n];
    c        =( -xs[n]*ys[n]/3.)/jac[n];
    d        =(a+xs[n]*xs[n]/3.)/jac[n];
    a        = a/jac[n];

    transport.getDViscosityDT   (1,&p,&t,&dmudT);
    transport.getDConductivityDT(1,&p,&t,&dkdT );

    f        = r*gm1/rGas;
    dTdq[0]  =-f*(e-(u*u+v*v));
    dTdq[1]  =-f*u;
    dTdq[2]  =-f*v;
    dTdq[3]  = f;
    dTdq[4]  = 0.;

    for (int k=0; k<nq; k++) dmudq[k] = dmudT*dTdq[k];

    chi      = rnu/mu;
    chi3     = chi*chi*chi;
    fv1      = chi3/(chi3+cv1*cv1*cv1);
    mut      = rnu*fv1;
    kt       = mut*rGas*ggm1/PrnT;
    fn       =(mu+rnu)/sigma;
    if (rnu < 0.){
      fn     =(cn1+chi3)/(cn1-chi3);
      fn     =(mu+rnu*fn)/sigma;
      mut    = 0.;
      kt     = 0.;
    }

    dfv1dchi  = 3.*chi*chi/(chi3+cv1*cv1*cv1)*(1.-fv1);

    for (int k=0; k<nq; k++) dchidq[k] =-chi*dmudq[k]/mu;
    dchidq[4] += 1./mu;
 
    for (int k=0; k<nq; k++) dfv1dq[k] = dfv1dchi*dchidq[k];

    for (int k=0; k<nq; k++) dmutdq[k] = rnu*dfv1dq[k];
    dmutdq[4] += fv1;

    dfndchi   = 0.;
    if (rnu < 0.) {
       dfndchi = 3.*chi*chi/(cn1-chi3)*(1.+fn);
       for (int k=0; k<nq; k++) dmutdq[k] = 0.;
    }

    for (int k=0; k<nq; k++) dfndq[k] = dfndchi*dchidq[k];

    for (int k=0; k<nq; k++) dwdq[k] = 1./sigma*(dmudq[k]+rnu*dfndq[k]);
    dwdq[4]  += fn/sigma;

    for (int k=0; k<nq; k++) dktdq[k] = rGas*ggm1/PrnT*dmutdq[k];
    
    mu       += mut;
    ku       += kt;
    for (int k=0; k<nq; k++) dmudq[k] += dmutdq[k];
    for (int k=0; k<nq; k++) dkdq[k]   = dkdT*dTdq[k]+dktdq[k];

    dumudq[0] = u*dmudq[0]-r*u*mu;
    dumudq[1] = u*dmudq[1]+r  *mu;
    dumudq[2] = u*dmudq[2]       ;
    dumudq[3] = u*dmudq[3]       ;
    dumudq[4] = u*dmudq[4]       ;

    dvmudq[0] = v*dmudq[0]-r*v*mu;
    dvmudq[1] = v*dmudq[1]       ;
    dvmudq[2] = v*dmudq[2]+r  *mu;
    dvmudq[3] = v*dmudq[3]       ;
    dvmudq[4] = v*dmudq[4]       ;

    A[iA   ] = 0.;
    A[iA+ 1] = 0.;
    A[iA+ 2] = 0.;
    A[iA+ 3] = 0.;
    A[iA+ 4] = 0.;

    A[iA+ 5] = dmudq[0]*(b*ui+c*vi);
    A[iA+ 6] = dmudq[1]*(b*ui+c*vi);
    A[iA+ 7] = dmudq[2]*(b*ui+c*vi);
    A[iA+ 8] = dmudq[3]*(b*ui+c*vi);
    A[iA+ 9] = dmudq[4]*(b*ui+c*vi);

    A[iA+10] = dmudq[0]*(c*ui+d*vi);
    A[iA+11] = dmudq[1]*(c*ui+d*vi);
    A[iA+12] = dmudq[2]*(c*ui+d*vi);
    A[iA+13] = dmudq[3]*(c*ui+d*vi);
    A[iA+14] = dmudq[4]*(c*ui+d*vi);

    A[iA+15] = dumudq[0]*(b*ui+c*vi)+dvmudq[0]*(c*ui+d*vi)+dkdq[0]*a*Ti;
    A[iA+16] = dumudq[1]*(b*ui+c*vi)+dvmudq[1]*(c*ui+d*vi)+dkdq[1]*a*Ti;
    A[iA+17] = dumudq[2]*(b*ui+c*vi)+dvmudq[2]*(c*ui+d*vi)+dkdq[2]*a*Ti;
    A[iA+18] = dumudq[3]*(b*ui+c*vi)+dvmudq[3]*(c*ui+d*vi)+dkdq[3]*a*Ti;
    A[iA+19] = dumudq[3]*(b*ui+c*vi)+dvmudq[4]*(c*ui+d*vi)+dkdq[4]*a*Ti;
    
    A[iA+20] = dwdq[0]*a*nui;
    A[iA+21] = dwdq[1]*a*nui;
    A[iA+22] = dwdq[2]*a*nui;
    A[iA+23] = dwdq[3]*a*nui; 
    A[iA+24] = dwdq[4]*a*nui;  
  }
}
