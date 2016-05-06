#ifndef included_ship
#define included_ship
class ship {
	public:
	
	float LifeVal, RangeVal;
	int MobilityVal, AttackVal;
	
	//Possible Moves of a Ship
	virtual void Move(float a) =0;
	//Possible Actions of a Ship
	virtual void Action(int i)=0;
	//Attack of a Ship
	virtual void Attack(float a)=0;
	//Deffence of a Ship
	virtual void Deffend(float a)=0;
	//Special Captin Abilities
	virtual void Captin(int i)=0;
	//Special Ship Abilities
	virtual void Shipa(int i)=0;
	
};
#endif
