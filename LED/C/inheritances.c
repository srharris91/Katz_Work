#include <iostream>
using namespace std;
class BigDude {
	private:
		int c;
	protected:
		int a,b;
	public:
		int e;
		void set_values(int g) {
			a=b=g;
			e=g;
			}
		friend class friendlydude;
	};
class LittleDude1: public BigDude {//inheritance1
	public:
		int area() { return (a*b);}
	};
class LittleDude2: public BigDude {//inerheritance2
	public:
		int area() {return(.5*a*b);}
	};
class friendlydude {
	public:
	int getfromBigDude (BigDude z) {	
		int h=z.a,i=z.b,j=z.c;
		return h;
	}
};
	int main() {
		LittleDude1 rec;
		LittleDude2 tri;
		rec.set_values(4);
		tri.set_values(4);
		cout<<rec.area()<<endl;
		cout<<tri.area()<<endl;
		cout<<tri.e<<" "<<endl;
		
		// for friend access my first try
		BigDude y;
		y.set_values(24);
		friendlydude w;
		cout<<w.getfromBigDude(y)<<endl;
		return 0;}
