// preprocessor test stuff.  Shorten compiler time.  #command
#ifndef __cplusplus
#error a c++ compiler is required you dumby
#endif

#include <iostream>
//~ #include <string>
using namespace std;//c++ code

#include "headerfileguy.h" //includes header file

#define g 5 //comment control this line  to control the ifdef functions 
#define str(x) #x
#define glue(x,y) x##y
#ifdef g
struct A {
	int h[g];
	string func(){return ("g is defined yo");}}; //c++ code
#endif

#ifndef g
#define g 10
struct A {
	int h[g];
	string func(){return ("g is not defined man");}}; //c++ code
#endif
//macro is named by _cplusplus this is the default macro named for all 
//c++ compilers
int main() {
	cout<<g<<endl;
	A C;
	cout<<C.func()<<endl;
	cout<<str(IDontBelieveThis)<<endl;
	glue(co,ut)<<"instead of writing cout"<<endl;
	int yu=29;
	glue(co,ut)<<glue(y,u)<<endl;
#line 1 "this is a new filename.c++"
//~ int a[]=; // comment control to see how #line works
	cout<<"This is the line number of this file "<<__LINE__<<endl;
	cout<<"This is the filename "<<__FILE__<<endl;
	cout<<"This is the data of compilation "<<__DATE__<<endl;
	cout<<"This is the time of compilation "<<__TIME__<<endl;
	cout<<"This is the cplusplus macro integer compiler number "
	<<__cplusplus<<endl;
	cout<<"CorC++ standard? "<<__STDC__<<endl; //if defined to 1, 
	//the implementation conforms to the C standard
	return 0;}
