//input and output stuff
#include <fstream>
#include <iostream>
using namespace std;
int main(){
	ofstream myfile;
	myfile.open ("outputexample.csv",ios::app);//append at end of file
	myfile<<"I'm going to write this to a file"<<endl;
	myfile<<"and maybe this too"<<endl;
	myfile.close();
	
	ofstream my2file ("output2example.csv",ios::app);
	string A="yo this is my second file";
	my2file<<A<<endl;
	
	ofstream bi ("example.bin");
	bi<<4<<5<<6<<7<<endl; //output for example.bin file
	
	ifstream my3file ("example.bin",ios::in | ios::binary | ios::ate);
	ifstream:: pos_type size; 
	char *memoryblockguy;
	size=my3file.tellg();//could have done=int size=int(my3file.tellg())
	//if the file is under 2 GB
	cout<<"how many bytes in file = "<<size<<endl;
	memoryblockguy = new char [size];//dynamic memory large enough
	// for the entire block
	my3file.seekg(0,ios::beg);
	my3file.read(memoryblockguy,size);//has to be read(char*,int)
	cout<<memoryblockguy[2]<<endl;
	delete[] memoryblockguy;
	cout<<"do this ? \? ffffffffffffffffffffffffffffffffffffff\r";
	cout<<"return this to the start"<<endl;
	cout<<"\n\n"<<endl;
	int O;
	cin>>O;
	my2file<<"\a \n\n"<<endl;
	my2file.close();
	my3file.close();
	return 0;}
