#include "Strand2dFCManager.h"


void Strand2dFCManager::mgMap(){
  int npos=0,level=0; bool down=true;
  Array1D<int> kcyc(nLevels); kcyc.set(0);
  while (true){
    //advance the solution in time on the current multigrid level
    if (down){
      if (level < nLevels-1){
	kcyc(level) += 1;
	level++;
	npos++;
      }
      else down = false;
    }
    else{
      if (level > 0){ //correct a fine level
	level--;
	npos++;
      }
      if (level == 0){
	down = true;
	break; //if back on fine level, cycle is done
      }
      if (kcyc(level) < mgCycle) down = true;
      else kcyc(level) = 0;
    }}

  mgLevel.allocate(npos+2);
  mgMode.allocate(npos+2);

  npos = 0; down = true; level = 0;
  mgLevel(npos) = level;
  mgMode(npos) = 1;
  kcyc.set(0);
  while (true){
    //advance the solution in time on the current multigrid level
    if (down){
      if (level < nLevels-1){
	kcyc(level) += 1;
	level++;
	npos++;
	mgLevel(npos) = level;
	mgMode(npos)  = 2;
      }
      else{
	down = false;
	if (nLevels > 1) mgMode(npos) = 3;
      }}
    else{
      if (level > 0){ //correct a fine level
	level--;
	npos++;
	mgLevel(npos) = level;
	mgMode(npos)  = 5;
      }
      if (level == 0){
	down = true;
	if (nLevels > 1) mgMode(npos) = 6;
	break; //if back on fine level, cycle is done
      }
      if (kcyc(level) < mgCycle){
	down = true;
	mgMode(npos) = 4;
      }
      else kcyc(level) = 0;
    }}

  mgLevel(npos+1) = 0;
  mgMode(npos+1) = 0;
  kcyc.deallocate();

  /*
  int n=0,mode=mgMode(0); level=mgLevel(0);
  do{
    cout << level << "\t" << mode << endl;
    n++;
    level = mgLevel(n);
    mode  = mgMode(n);
  } while(mode != 0);
  cout << level << "\t" << mode << " exit." << endl;
  exit(0);
  */
}
