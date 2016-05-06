#include "Ship.h"
#include <iostream>
using namespace std;
class Sloop: public ship {
	public:
	Sloop(){
		LifeVal=100.;
		RangeVal=4.;
		MobilityVal=3;
		AttackVal=8;
	};
	
	void Move(float a) {
		
		cout << "Move S";
		
	};
	
	void Action(int i) {
		cout << "Action S";
		cout << LifeVal << '\n';
	};
	void Attack(float a) {
		cout << "Attack S";
	};
	void Deffend(float a) {
		cout << "Deffend S";
	};
	void Captin(int i) {
		cout << "Cap S";
	};
	void Shipa(int i) {
		cout << "Ship S";
	};
};
