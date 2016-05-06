// LED_FirstDWave simple
//~ Note:  Function is of the form du/dt+a*du/dx=f(x)
//~  this takes a reconstructed with shockwave form giving a 2nd order accurate values
#include <iostream>//normal other stuff
#include <fstream> //input from file and output to file
#include <math.h>  //sin and cos
#include <string>  //string variables
//~ #include <stdlib.h>//string to double converter atof function
using namespace std;

// user input source and initial value functions
double source(double x) {//source function
	//~ return 0;};//no source term
	return (-1./5.*(x-5.));}; //input function for source term
double IV(double x) {//initial value function
	return(cos(2.*x));};
int main() {
	ofstream output,test;
	ifstream input;
	// input data and read from file
	output.open ("LEDOutput.csv");
	input.open ("LEDInput.csv");
	string A; //dummy string variable
	double l,dx,time,dt,a;
	input>>A>>A>>A>>A>>A;
	cout <<A<<endl;
	input >> l>>dx>>time>>dt>>a;
	cout << "length = "<<l << "  dx = "<<dx << "  time = "<<time
		<< "  dt = "<<dt <<"  a = "<< a << endl;
	// calc using LED method
	int dxi=int(l/dx+1),dti=int(time/dt+1);//spatial and time nodes
	cout<<"spatial nodes = "<<dxi<<"  time nodes = "<<dti<<endl;
	double u[dti][dxi+3];//initialize u array with ghost nodes, 2 on each side of bar
	double urr,url,ull,ulr,dphalf,dmhalf,umhalf,uphalf;//other LED variables
	double xmh,xph,trapS;
	double sofa=(a>-1)?1:-1;//sign of a
	cout<<sofa<<endl;
	for (int i=0;i<dti;i++){//time
		for (int j=0;j<2;j++){//2 left spatial ghost nodes
			u[i][j]=double (0.0);}
		for (int j=dxi+2;j<(dxi+4);j++){//2 right spatial ghost nodes
			u[i][j]=double(0.0);}
			}
	for (int j=2;j<dxi+2;j++){//bar spatial nodes initial values
		u[0][j]=double(IV(double(j-2)*dx));}
	for (int i=0;i<(dti-1);i++) {//time jump
		for (int j=2;j<(dxi+2);j++) {//spatial jump
			urr=u[i][j+1]-sofa*dx/2.*(u[i][j+2]-u[i][j])/(2.*dx);
			url=u[i][j]+sofa*dx/2.*(u[i][j+1]-u[i][j-1])/(2.*dx);
			ull=u[i][j-1]+sofa*dx/2.*(u[i][j-2]-u[i][j])/(2.*dx);
			ulr=u[i][j]-sofa*dx/2.*(u[i][j+1]-u[i][j-1])/(2.*dx);
	//~ no absolute value in dphalf or dmhalf to handle negative a
			dphalf=.5*(a)*(urr-url);
			dmhalf=.5*(a)*(ulr-ull);
			umhalf=.5*fabs(a)*(u[i][j]+u[i][j-1])-dmhalf;
			uphalf=.5*fabs(a)*(u[i][j]+u[i][j+1])-dphalf;
	//~ single trapezoidal of source term from x-1/2 to x+1/2
			xmh=(double (j-2)-.5)*dx;
			xph=(double (j-2)+.5)*dx;
			trapS=.5*(xph-xmh)*(source(xph)+source(xmh));
			u[i+1][j]=-a*dt*(uphalf-umhalf)/dx+u[i][j]+trapS*dt/dx;
		}
	}
	output.precision(17); //output 17 digits of precision
	for (int i=0;i<(dti);i++){
		for (int j=2;j<(dxi+2);j++){
			output<<double(u[i][j])<<" ";
		}
		output<<endl;
	}
	input.close();
	output.close();
return 0;
}
