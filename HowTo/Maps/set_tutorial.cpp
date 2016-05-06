#include <iostream>
#include <set>

using namespace std;

int main(){
    set<char> abc;      // initialize the emtpy set
    abc.insert('b');    // insert values
    abc.insert('c');    // insert values
    abc.insert('a');    // insert values
    abc.insert('a');    // insert values

//    set<char>::iterator it; // initialize iterator for loop
    for (auto&& it=abc.begin(); it!=abc.end(); ++it){  // loop through all values
        cout<<" "<<*it<<endl;}                  // output the values contained in set


    set<int> a123;
    a123.insert(2);
    a123.insert(1);
    a123.insert(2);
    a123.insert(3);
    a123.insert(4);
//    set<int>::iterator it2; // initialize iterator for loop
    for (auto it=a123.begin(); it!=a123.end(); ++it){  // loop through all values
        cout<<" "<<*it<<endl;}                  // output the values contained in set


return 0;}
