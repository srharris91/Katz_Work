// bad alloc trial exceptions and stuff
#include <iostream>
//~ #include <exception>
using namespace std;


int main()
{
	int * a=new int[10000000000];
	cin>>a[5];
	return 0;
}
