#include "solutionPoints1D.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
using namespace std;


// returns solution points based on desired order and spacing in 1D

void solutionPoints1D(const int& order,
		      const int& spacing,
		      double* r)
{
  if      (spacing == 0){ //equal spacing
    int m=0;
    if      (order == 0){
      r[0] = 0.;
      }
    else if (order == 1){
      r[0] =-1.;
      r[1] = 1.;
    }
    else if (order == 2){
      r[0] =-1.;
      r[1] = 1.;
      r[2] = 0.;
    }
    else if (order == 3){
      r[0] =-1.;
      r[1] = 1.;
      r[2] =-1./3.;
      r[3] = 1./3.;
    }
    else if (order == 4){
      r[0] =-1.;
      r[1] = 1.;
      r[2] =-.5;
      r[3] = 0.;
      r[4] = .5;
    }
    else{
      cout << "\nChoose a lower order in solutionPoints.C\n";
      exit(0);
    }
  }


  // else an error
  else{
    cout << "\nPoint spacing not recognized in solutionPoints1D.C\n";
    exit(0);
  }
}
