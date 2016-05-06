#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::fdCoeffInviscid()
{
  // inviscid coefficients of various orders
  if (strandOrder == 1 || strandOrder == 2){ //D_1^(2,1,2)
    // dimensions
    int nBnd=2,nMax=3;
    if (nStrandNode < 2*nBnd){
      cout << "\nPlease choose nStrandNode greater than "
	   << 2*nBnd << "." << endl;
      exit(0);
    }

    icn2 = new int*[nStrandNode];
    for (int j=0; j<nStrandNode; j++) icn2[j] = new int[2];

    icn2[0][0] = 0; //root stencils
    icn2[0][1] = 2;

    icn2[1][0] = 0;
    icn2[1][1] = 3;

    for (int j=nBnd; j<nStrandNode-nBnd; j++){ //interior stencils
      icn2[j][0] = j-1;
      icn2[j][1] = 3;
    }

    int i=nBnd-1,j1,nj,l;
    for (int j=nStrandNode-nBnd; j<nStrandNode; j++){ //tip stencils
      j1         = icn2[i  ][0];
      nj         = icn2[i--][1];
      icn2[j][0] = nStrandNode-(j1+nj);
      icn2[j][1] = nj;
    }

    // coefficients
    icn1 = new double*[nStrandNode];
    for (int j=0; j<nStrandNode; j++) icn1[j] = new double[icn2[j][1]];

    icn1[0][0] =-1.; //root coefficients
    icn1[0][1] = 1.;

    icn1[1][0] =-.5;
    icn1[1][1] = 0.;
    icn1[1][2] = .5;

    for (int j=nBnd; j<nStrandNode-nBnd; j++){ //interior coefficients
      icn1[j][0] =-.5;
      icn1[j][1] = 0.;
      icn1[j][2] = .5;
    }

    i = nBnd-1;
    for (int j=nStrandNode-nBnd; j<nStrandNode; j++){ //tip coefficients
      j1         = icn2[i][0];
      nj         = icn2[i][1];
      l          = nj-1;
      for (int m=0; m<nj; m++) icn1[j][m] =-icn1[i][l--];
      i--;
    }

    // inverse of norm matrix (first element) needed for penalty BC term
    Pinv0 = 2.;
  }


  else if (strandOrder == 3){ //D_1^(4,2,3)
    // dimensions
    int nBnd=4,nMax=6;
    if (nStrandNode < 2*nBnd){
      cout << "\nPlease choose nStrandNode greater than "
	   << 2*nBnd << "." << endl;
      exit(0);
    }
    icn2 = new int*[nStrandNode];
    for (int j=0; j<nStrandNode; j++) icn2[j] = new int[2];

    icn2[0][0] = 0; //root stencils
    icn2[0][1] = 4;

    icn2[1][0] = 0;
    icn2[1][1] = 3;

    icn2[2][0] = 0;
    icn2[2][1] = 5;

    icn2[3][0] = 0;
    icn2[3][1] = 6;

    for (int j=nBnd; j<nStrandNode-nBnd; j++){ //interior stencils
      icn2[j][0] = j-2;
      icn2[j][1] = 5;
    }

    int i=nBnd-1,j1,nj,l;
    for (int j=nStrandNode-nBnd; j<nStrandNode; j++){ //tip stencils
      j1         = icn2[i  ][0];
      nj         = icn2[i--][1];
      icn2[j][0] = nStrandNode-(j1+nj);
      icn2[j][1] = nj;
    }

    // coefficients
    icn1 = new double*[nStrandNode];
    for (int j=0; j<nStrandNode; j++) icn1[j] = new double[icn2[j][1]];

    icn1[0][0] =-24.0/17.0; //root coefficients
    icn1[0][1] = 59.0/34.0;
    icn1[0][2] =-4.0/17.0;
    icn1[0][3] =-3.0/34.0;

    icn1[1][0] =-.5;
    icn1[1][1] = 0.;
    icn1[1][2] = .5;

    icn1[2][0] = 4.0/43.0;
    icn1[2][1] =-59.0/86.0;
    icn1[2][2] = 0.;
    icn1[2][3] = 59.0/86.0;
    icn1[2][4] =-4.0/43.0;

    icn1[3][0] = 3.0/98.0;
    icn1[3][1] = 0.0;
    icn1[3][2] =-59.0/98.0;
    icn1[3][3] = 0.0;
    icn1[3][4] = 32.0/49.0;
    icn1[3][5] =-4.0/49.0;

    for (int j=nBnd; j<nStrandNode-nBnd; j++){ //interior coefficients
      icn1[j][0] = 1./12.;
      icn1[j][1] =-2./3.;
      icn1[j][2] = 0.;
      icn1[j][3] = 2./3.;
      icn1[j][4] =-1./12.;
    }

    i = nBnd-1;
    for (int j=nStrandNode-nBnd; j<nStrandNode; j++){ //tip coefficients
      j1         = icn2[i][0];
      nj         = icn2[i][1];
      l          = nj-1;
      for (int m=0; m<nj; m++) icn1[j][m] =-icn1[i][l--];
      i--;
    }

    // inverse of norm matrix (first element) needed for penalty BC term
    Pinv0 = 48./17.;
  }


  else if (strandOrder == 4){ //D_1^(6,3,4)
    // dimensions
    int nBnd=6,nMax=9;
    if (nStrandNode < 2*nBnd){
      cout << "\nPlease choose nStrandNode greater than "
	   << 2*nBnd << "." << endl;
      exit(0);
    }
    icn2 = new int*[nStrandNode];
    for (int j=0; j<nStrandNode; j++) icn2[j] = new int[2];

    icn2[0][0] = 0; //root stencils
    icn2[0][1] = 6;

    icn2[1][0] = 0;
    icn2[1][1] = 6;

    icn2[2][0] = 0;
    icn2[2][1] = 6;

    icn2[3][0] = 0;
    icn2[3][1] = 7;

    icn2[4][0] = 0;
    icn2[4][1] = 8;

    icn2[5][0] = 0;
    icn2[5][1] = 9;

    for (int j=nBnd; j<nStrandNode-nBnd; j++){ //interior stencils
      icn2[j][0] = j-3;
      icn2[j][1] = 7;
    }

    int i=nBnd-1,j1,nj,l;
    for (int j=nStrandNode-nBnd; j<nStrandNode; j++){ //tip stencils
      j1         = icn2[i  ][0];
      nj         = icn2[i--][1];
      icn2[j][0] = nStrandNode-(j1+nj);
      icn2[j][1] = nj;
    }

    // coefficients
    icn1 = new double*[nStrandNode];
    for (int j=0; j<nStrandNode; j++) icn1[j] = new double[icn2[j][1]];

    //from Zingg or Hiener?
    icn1[0][0] =-1.58253351893912; //root coefficients
    icn1[0][1] = 2.01907954184678;
    icn1[0][2] =-0.0843163113292793;
    icn1[0][3] =-0.536193127701663;
    icn1[0][4] = 0.161684616699636;
    icn1[0][5] = 0.0222787994236452;

    icn1[1][0] =-0.458809900385693;
    icn1[1][1] = 0.;
    icn1[1][2] = 0.254765670523599;
    icn1[1][3] = 0.323801992286135;
    icn1[1][4] =-0.117851494214601;
    icn1[1][5] =-0.00190626820943978;

    icn1[2][0] = 0.0424505102668142;
    icn1[2][1] =-0.564459609000369;
    icn1[2][2] = 0.;
    icn1[2][3] = 0.462252551334071;
    icn1[2][4] = 0.121080781999262;
    icn1[2][5] =-0.0613242345997787;

    icn1[3][0] = 0.136564657585370;
    icn1[3][1] =-0.362925297008148;
    icn1[3][2] =-0.233843378739815;
    icn1[3][3] = 0.;
    icn1[3][4] = 0.483843378739815;
    icn1[3][5] =-0.0370747029918517;
    icn1[3][6] = 0.0134353424146296;

    icn1[4][0] =-0.0560323304134400;
    icn1[4][1] = 0.179732131522153;
    icn1[4][2] =-0.0833439126571030;
    icn1[4][3] =-0.658351318183742;
    icn1[4][4] = 0.;
    icn1[4][5] = 0.764244001523423;
    icn1[4][6] =-0.164529643265203;
    icn1[4][7] = 0.0182810714739114;

    icn1[5][0] =-0.00694238335502234;
    icn1[5][1] = 0.00261409556859433;
    icn1[5][2] = 0.0379557544348303;
    icn1[5][3] = 0.0453604560017656;
    icn1[5][4] =-0.687193214766786;
    icn1[5][5] = 0.;
    icn1[5][6] = 0.739709139060752;
    icn1[5][7] =-0.147941827812150;
    icn1[5][8] = 0.0164379808680167;

    /*
    // from Mattson
    icn1[0][0] =-21600./13649.; //root coefficients
    icn1[0][1] = 104009./54596.;
    icn1[0][2] = 30443./81894.;
    icn1[0][3] =-33311./27298.;
    icn1[0][4] = 16863./27298.;
    icn1[0][5] =-15025./163788.;

    icn1[1][0] =-104009./240260.;
    icn1[1][1] = 0.;
    icn1[1][2] =-311./72078.;
    icn1[1][3] = 20229./24026.;
    icn1[1][4] =-24337./48052.;
    icn1[1][5] = 36661./360390.;

    icn1[2][0] =-30443./162660.;
    icn1[2][1] = 311./32532.;
    icn1[2][2] = 0.;
    icn1[2][3] =-11155./16266.;
    icn1[2][4] = 41287./32532.;
    icn1[2][5] =-21999./54220.;

    icn1[3][0] = 33311./107180.;
    icn1[3][1] =-20229./21436.;
    icn1[3][2] = 485./1398.;
    icn1[3][3] = 0.;
    icn1[3][4] = 4147./21436.;
    icn1[3][5] = 25427./321540.;
    icn1[3][6] = 72./5359.;

    icn1[4][0] =-16863./78770.;
    icn1[4][1] = 24337./31508.;
    icn1[4][2] =-41287./47262.;
    icn1[4][3] =-4147./15754.;
    icn1[4][4] = 0.;
    icn1[4][5] = 342523./472620.;
    icn1[4][6] =-1296./7877.;
    icn1[4][7] = 144./7877.;

    icn1[5][0] = 15025./525612.;
    icn1[5][1] =-36661./262806.;
    icn1[5][2] = 21999./87602.;
    icn1[5][3] =-25427./262806.;
    icn1[5][4] =-342523./525612.;
    icn1[5][5] = 0.;
    icn1[5][6] = 32400./43801.;
    icn1[5][7] =-6480./43801.;
    icn1[5][8] = 720./43801.;
    */

    for (int j=nBnd; j<nStrandNode-nBnd; j++){ //interior coefficients
      icn1[j][0] =-1./60.;
      icn1[j][1] = 3./20.;
      icn1[j][2] =-3./4.;
      icn1[j][3] = 0.;
      icn1[j][4] = 3./4.;
      icn1[j][5] =-3./20.;
      icn1[j][6] = 1./60.;
    }

    i = nBnd-1;
    for (int j=nStrandNode-nBnd; j<nStrandNode; j++){ //tip coefficients
      j1         = icn2[i][0];
      nj         = icn2[i][1];
      l          = nj-1;
      for (int m=0; m<nj; m++) icn1[j][m] =-icn1[i][l--];
      i--;
    }

    // inverse of norm matrix (first element) needed for penalty BC term
    Pinv0 = 43200./13649.;
  }

  else{
    cout << "\nThis strand order not yet available for inviscid terms." << endl;
    exit(0);
  }

  /*
  int nj,j1;
  double sum;
  for (int j=0; j<nStrandNode; j++){
    j1 = icn2[j][0];
    nj = icn2[j][1];
    cout << "\nNode: " << j << " " << endl;
    sum = 0.;
    for (int m=0; m<nj; m++){
      cout << j1++ << " " << icn1[j][m] << endl;
      sum += icn1[j][m];
    }
    cout << sum << endl;
  }
  exit(0);
  */
}
