// Computer 2 copy of fortran code to c++
#include <iostream>
#include <math.h>
using namespace std;
double modsec(double& f,double d,double epd,double Re) {
	while (1==1) {
	f=f-d*f*(1./sqrt(f)+2.*log10(epd/3.7+2.51/(Re*sqrt(f)))) \
			/(1./sqrt(f+d*f)+2.*log10(epd/3.7+2.51/(Re*sqrt(f+d*f)))- \
			(1./sqrt(f)+2.*log10(epd/3.7+2.51/(Re*sqrt(f)))));
	if (fabs(1./sqrt(f)+2.*log10(epd/3.7+2.51/(Re*sqrt(f))))<=0.000001) {break;}
	cout<<f<<endl; //so the user can see the value of each iteration;
	return(f);};
}
	int main() {
		double f,d=.01,epd,Re;
		cout<<"please input f initial guess, epd, and Re"<<endl;
		cin>>f>>epd>>Re;
		modsec(f,d,epd,Re);
		cout<<f<<endl;
		return 0;}
