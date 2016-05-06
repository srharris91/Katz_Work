#include <iostream>
using namespace std;
int main () {
	int n;
	cout<<"Enter value for n : ";
	cin>>n;
	if (n>0) {
		while (n>0)
			{cout << n<<endl;n--;}
			cout<<"Fire!"<<endl;}
	else if (n<0) {
		while (n<0) {
			cout<<n<<endl;n++;}
			cout<<"Fire!"<<endl;}
	else
		{cout<<"n=0 \n";}
		
	do {
		cout<<"yo"<<++n<<endl;}
		while (n<10);
	int x;
	cout<<"input x value (1,2,or 3)"<<endl;
	cin>>x;
	switch (x) {
		case 1:
			cout<<"x= "<<x<<endl;
			break;
		case 2:
			cout<<"x= "<<x<<endl;
			break;
		case 3:
			cout<<"x= "<<x<<endl;
			break;
		default:
			cout<<"x!=1,2,or 3"<<endl;
			cout<<"x= "<<x<<endl;
			break;}
	
return 0; }

