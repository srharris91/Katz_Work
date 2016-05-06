// play with dynamic memory and stuff
#include <iostream>
using namespace std;

struct structthingy {
	int y,z,b,c;
	structthingy *a;
	int functionthingy() {
	a=this;
}
};

int main() {
	int *dynamicarray;
	int n;
	cout<<"How large of an array?"<<endl;
	cin>>n;
	dynamicarray=new int [n];
	cout<<"Please input array"<<endl;
	for (int i=0;i<n;i++) {
		cin>>dynamicarray[i];}
	for (int i=0;i<n;i++) {
		cout<<dynamicarray[i];}
		//~ structthingy.functionthingy();
		
	return 0; }
