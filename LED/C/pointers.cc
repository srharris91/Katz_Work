// pointers
#include <iostream>
using namespace std;

void increase(void*data,int psize)
{
	if (psize == sizeof(char))
	{
		char*pchar=(char*)data; ++*pchar;
	}
	else if (psize == sizeof(int))
	{
		int*pint=(int*)data; ++(*pint);
	}
}

struct A {
	int j,k,l;
	A* m(){return(this);};};

int main()
{
	char a='x';
	int b=11;
	int*c=&b;
	increase(&a,sizeof(a));
	increase(&b,sizeof(b));
	cout<<a<<","<<b<<endl;
	int arraystuff[]={1,2,3,4,5,6,7,8,9};
	double secondarray[]={1.,2.23456,3.45,4.5,5.,6.};
	cout<<sizeof(arraystuff)/sizeof(int)<<endl;
	cout<<sizeof(secondarray)/sizeof(double)<<sizeof(double)<<endl;
	cout<<c<<" "<<&b<<endl;
	cout<<*c<<" "<<b<<endl;
	A B;
	B.j=4;
	cout<<B.j<<endl;
	A *C,*D;
	C=&B;
	cout<<B.m()<<" "<<C<<endl;
	// * is value pointed by
	// & is address of
	return 0;
}
