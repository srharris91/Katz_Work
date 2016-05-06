//type casting stuff
#include <iostream>
using namespace std;

class Basedude { virtual void functionname() {}};
class Deriveddude: public Basedude {int k;};


int main () {
	short a=32767;// 2 bytes... maximum is 32767 minimum is -32678
	int b; // 4 bytes
	double c; // 8 bytes
	b=a; // convert implicitly
	c=double(a); // function look to convert explicitly
	cout<<a<<" short "<<b<<" int "<<c<<" double"<<endl;
	
	double *w;
	int *v=&b;
	w=(double*) v;// must explicity convert pointer types
	cout<<w<<" "<<v<<endl;
	
	Deriveddude z,*Z;
	Basedude X,*x=dynamic_cast<Basedude*>(&z);// dynamic cast
	cout<<&z<<" "<<x<<endl;
	Z=static_cast<Deriveddude*>(&X);//static cast
	cout<<&X<<" "<<Z<<endl;
	
	
	return 0;}
