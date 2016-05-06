// my first program in C++
// designates a comment line...
/* this does too
  and it can extend to many
  many
  lines */
#include <iostream>
#define di .1234567890123456789
using namespace std;
int main ()
{
	double b=3.123456789012345678901234567890l,a,c,d;
	int x=2;
	string f="really long string of stuff";
	a=c=2.0;
	cout << a << b << c <<endl;
	b=5.+(d=10.0);
	cout << b << d << endl;
	//~ bool a=false;
	//~ cout << "Hello World!  \n";
	//~ cout << "yo\n";
	//~ cout << a << endl;
	cout.precision(15); // show double precision length
	cout << 3.123456789012345678901234567890L << endl; // Long Double precision
	cout << 3.123456789012345678901234567890f << endl; // floating point single precision
	cout << 0113 << endl;
	cout << 0x4b << endl;
	cout << di << endl;
	cout << "beep yo \a" << endl;
	// see *= does, can be done with += or -= or /=
	cout << x << endl;
	x *= 2;
	cout << x << endl;
	cout << f << endl;
	cin>>d;
	cout<<d<<endl;
	return 0;
}
