#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
using namespace std;
#include "headerfileguy.h"
int main() {
	int y=5,g=2,yoyo=10;
	stringstream a;
	//~ a.str(""); //clears a.str() value
	a.clear();
	a<<"this is awesome"<<y<<g<<yoyo<<endl; //input files into a
	cout<<a.str()<<endl; //show output
	string b=a.str(); //copy output into a string b
	a.clear(); //clear flags
	a.str("");
	cout<<a.str()<<endl;
	cout<<"well that was a waste"<<endl;
	cout<<b<<endl;// output copy string
	ifstream plotMesh("stringinputthingy.txt");

	string stry3;
	getline(plotMesh,stry3);
	
	cout<<stry3<<endl;
	cout<<Getmin(5,4)<<"\n"<<endl;
	
	//built in macros
	cout<<__LINE__<<endl; 
	cout<<__FILE__<<endl;
	cout<<__DATE__<<endl;
	cout<<__TIME__<<endl;
	cout<<__cplusplus<<endl;
	
	
	plotMesh.close();
	return 0;}
