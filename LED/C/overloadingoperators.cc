//structures classes and overloading operators
#include <iostream>
using namespace std;

struct cvec {
	int x,y,z;
	cvec(){x=0;y=0;z=0;};
	cvec(int a, int b, int c) {x=a;y=b;z=c;};
	cvec operator + (cvec param) {
	cvec temp;
		temp.x=x+param.x;
		temp.y=y+param.y;
		temp.z=z+param.z;
		return (temp);
	}
};
//~ cvec cvec::operator + (cvec param) {
	//~ 
	//~ cvec temp;
		//~ temp.x=x+param.x;
		//~ temp.y=y+param.y;
		//~ temp.z=z+param.z;
		//~ return (temp);
	//~ }
int main(){
	cvec a (3,1,2);
	cvec b (1,2,10);
	cvec c,d;
	c=a+b;
	d=a.operator+(b);
	cout<<c.x<<" "<<c.y<<" "<<c.z<<endl;
	cout<<d.x<<" "<<d.y<<" "<<d.z<<endl;
	return 0;
	}
