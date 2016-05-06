//polymorphism and stuff
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
		virtual int area(void)=0;//pure virtual function
		void printstuff (void){
			cout<<this->area()<<endl;}
	};
class LittleDude1: public BigDude {//inheritance1
	public:
		int area(void) { return (a*b);}
	};
class LittleDude2: public BigDude {//inerheritance2
	public:
		int area(void) {return(.5*a*b);}
	};

	int main() {
		LittleDude1 rec;
		LittleDude2 tri;
		BigDude *pointer1 = &rec;
		BigDude *pointer2 = &tri;
		pointer1->set_values(4);
		pointer2->set_values(4);
		cout<<rec.area()<<" "<<pointer1->area()<<endl;
		pointer1->printstuff();
		pointer2->printstuff();
		cout<<tri.area()<<" "<<pointer2->area()<<endl;
		cout<<tri.e<<" "<<pointer2->e<<endl;
		return 0;}
