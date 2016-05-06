#include <iostream>
#include <map>
#include <string>

using namespace std;

int main(){
    map<string,string> yo;
    yo["John"]="A+";
    cout<<yo["John"]<<endl;
    yo["John"]="A-";
    cout<<yo["John"]<<endl;



return 0;}
