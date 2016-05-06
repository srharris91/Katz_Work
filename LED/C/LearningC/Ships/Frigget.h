#include "Ship.h"
#include <iostream>
using namespace std;

class Frigget: public ship {
	public:
	Frigget(){
		LifeVal=100.;
		RangeVal=4.;
		MobilityVal=3;
		AttackVal=8;
	};
	
	void Move(float a) {
		
		cout << "Move";
		
	};
	
	void Action(int i) {
		cout << "Action";
	};
	void Attack(float a) {
		cout << "Attack";
	};
	void Deffend(float a) {
		cout << "Deffend";
	};
	void Captin(int i) {
		cout << "Cap";
	};
	void Shipa(int i) {
		cout << "Ship";
	};
};
