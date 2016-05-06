#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::limit(const int& j)
{
/* this calculates the limiter value for Strand2dFC
 * written by Shaun Harris at Utah State University
 * email: shaun.r.harris@gmail.com
 *  Sept. 16, 2014
 * 
 * 
 * 
 */

// standard limiter=0 from input.namelist
	if (limiter == 0) {
		for(int n=0; n<nSurfNode; n++){
			for (int k=0; k<nq; k++) {
				lim(n,j,k) = 1.;
			}
		}
	}
	
// limiter to be used limiter=1 from input.namelist
	else if (limiter == 1){ //calculate using gradient
		double a;
		double b;
		int n0;
		int n1;
		double limE;
		//initially set all lim values of j to 1
		for(int n=0; n<nSurfNode; n++){
			for (int k=0; k<nq; k++) {
				lim(n,j,k) = 1.;}}
		//find limiter value at each surfEdge
		for(int n=0; n<nSurfEdge; n++){
		    for (int k=0; k<nq; k++) {
				n0=surfEdge(n,0);
				n1=surfEdge(n,1);
				a=deltaS*qx(n1,j,k,0); 
				b=deltaS*qx(n0,j,k,0);	
				//~ a = q(surfEdge(n,3),j,k)-q(n1,j,k); // original a and b, can use this instead of a and b above
				//~ b = q(n0,j,k)-q(surfEdge(n,2),j,k);
				limE = 1.-((fabs((a-b)/max((fabs(a)+fabs(b)),dlim(k))))*
				(fabs((a-b)/max((fabs(a)+fabs(b)),dlim(k))))*
				(fabs((a-b)/max((fabs(a)+fabs(b)),dlim(k)))));
				//output warning and exit if bad limiter
				if (limE>1 | limE<0){
				    cout<<"WARNING: Bad Limiter value"<<endl;
				    cout<<"limE("<<n<<","<<j<<","<<k<<") = "<<limE<<" "<<"limiter "<<limiter<<endl;
				    exit(0);
					}
				//convert limiter value from Edge to Node using the minimum of 2 neighboring edges of each node
				lim(n0,j,k)=min(lim(n0,j,k),limE);
				lim(n1,j,k)=min(lim(n1,j,k),limE);
				//~ cout<<"lim("<<n0<<","<<j<<","<<k<<")= "<<lim(n0,j,k)<<endl;
				//~ cout<<"lim("<<n1<<","<<j<<","<<k<<")= "<<lim(n1,j,k)<<endl<<endl;
			}
		}
	}
	else {
		cout << "WARNING: input value for limiter is incorrect format  Limiter = "<<limiter<<endl;
		exit(0);
	}
	
}
