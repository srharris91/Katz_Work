// Class Ship
#include "Ships/shiplib.h"
#include <iostream>
using namespace std;

int main () {
	int nShips;
	cout << "Learning C .1"<<endl;
	cout << "How many ships?" <<endl;
	cin >> nShips;
    ship* S1 = new Frigget;
	ship* S2 = new Sloop;
	ship ** Frends = new ship*[nShips];
	for (int i=1; i<=nShips; i++) {
		int St;
		cout <<"Ship type 1 or 2?";
		cin >> St;
		switch(St) {
			case 1:
				cout <<"Frigget Created"<<endl;
				Frends[i] = S1;
				break;
			case 2:
				cout <<"Sloop Created"<<endl;
				Frends[i] = S2;
				break;
		}
	} 
	return 0;
}
