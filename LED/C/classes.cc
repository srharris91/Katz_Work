// classes and stuff
#include <iostream>
using namespace std;

class CTry {
	double x,y,z;
	public:
		//~ double setvalues (double a,double b,double c)
		//~ {
			//~ x=a;y=b;z=c;
		//~ }
		CTry(double,double,double);
		double volume() {return x*y*z;}
		int randomness(int a, int b, int c)
		{
			cout<<a<<" "<<x<<endl;
			return 0;
		}
		double a,b,c,d,e,f,g,h;
		struct yo {
			int i,j,k;
			double l,m;} n,o,p;
};
CTry::CTry (double a, double b, double c) { // Constructor of class CTry
	x=a;y=b;z=c; } // this only works when objects of this class are created
int main()
{
	// must input variables like this to have the constructor work correctly
	CTry cube1(1,2,3),cube2(2,3,4); 
	//~ cube1.setvalues(1,2,3);
	//~ cube2.setvalues(2,3,4);
	cout<<cube1.volume()<<endl<<cube2.volume()<<endl;
	cube1.randomness(9,8,7);
	cube2.randomness(9,8,7);
	cube1.a=20;
	cout<<cube1.a<<endl;
	cube1.b=5.45;
	cout<<cube1.b<<endl;
	cube1.n.i=98;
	cout<<cube1.n.i<<endl;
	return 0;
}
