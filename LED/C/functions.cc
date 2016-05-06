#include <iostream>
using namespace std;

double power(double a,int b) { //first function
double r=a;
// first value raised to the power of the second value
while (b>1) {
r=r*a;
b--;}
return(r);}

double factorial(double a) {
if (a>1)
return(a*factorial(a-1));
else
return(1);}
void writestatement () { // second function using void type
cout << "Write somethin' to the screen" << endl;
}
void duplicate (double& a, double& b, double& c) { 
//reference variables instead of the values using &
a*=2;b*=2;c*=2; }

//used to make functions code body written after int main ()
void preorder(int x,int y);

int main() {
// first value raised to the power of the second value
cout<<power(-4.,5)<<endl; 
writestatement();

double d=4, e=3, f=2;
duplicate (d,e,f);
cout<<d<<endl<<e<<endl<<f<<endl;
cout<<factorial(4)<<endl;
while(d<10)
	{
	d++;
	for(int i=4,i<10,i++) 
	{
	cout<<i<<endl;
	}
}
preorder(d,98);
return (0);}

//used to make functions code body written after int main ()
void preorder(int x,int y) { 
cout<< "just print values to screen "<<x<<" "<<y<<endl;}
