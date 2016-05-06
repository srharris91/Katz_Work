#include "Strand2dFCBlockSolver.h"


void Strand2dFCBlockSolver::fdCoeffDissipation()
{
  if (strandOrder == 1){ // Dissipation_2_1
    // dimensions
    int nBnd=3,nMax=2;
    if (nStrandNode < 2*nBnd){
      cout << "\nPlease choose nStrandNode greater than "
	   << 2*nBnd << "." << endl;
      exit(0);
    }
    
    dcn2 = new int*[nStrandNode+1];
    for (int j=0; j<nStrandNode+1; j++) dcn2[j] = new int[2];
    
    dcn2[0][0] = 0; //root stencils
    dcn2[0][1] = 2;
    
    dcn2[1][0] = 0;
    dcn2[1][1] = 2;
    
    dcn2[2][0] = 1;
    dcn2[2][1] = 2;
    
    for (int j=nBnd; j<nStrandNode+1-nBnd; j++){ //interior stencils
      dcn2[j][0] = j-1;
      dcn2[j][1] = 2;
    }
    
    int i=nBnd-1,j1,nj,l;
    for (int j=nStrandNode+1-nBnd; j<nStrandNode+1; j++){ //tip stencils
      j1         = dcn2[i  ][0];
      nj         = dcn2[i--][1];
      dcn2[j][0] = nStrandNode-(j1+nj);
      dcn2[j][1] = nj;
    }
    
    // coefficients
    dcn1 = new double*[nStrandNode+1];
    for (int j=0; j<nStrandNode+1; j++) dcn1[j] = new double[dcn2[j][1]];
    
    dcn1[0][0] = 1.; //root coefficients
    dcn1[0][1] =-1.;
    
    dcn1[1][0] =-1.;
    dcn1[1][1] = 1.;
    
    dcn1[2][0] =-1.;
    dcn1[2][1] = 1.;
    
    for (int j=nBnd; j<nStrandNode+1-nBnd; j++){ //interior coefficients
      dcn1[j][0] =-1.;
      dcn1[j][1] = 1.;
    }
    
    i = nBnd-1;
    for (int j=nStrandNode+1-nBnd; j<nStrandNode+1; j++){ //tip coefficients
      j1         = dcn2[i][0];
      nj         = dcn2[i][1];
      l          = nj-1;
      for (int m=0; m<nj; m++) dcn1[j][m] =-dcn1[i][l--];
      i--;
    }
  }


  else if (strandOrder == 2){ // Dissipation_4_2
    // dimensions
    int nBnd=5,nMax=6;
    if (nStrandNode < 2*nBnd){
      cout << "\nPlease choose nStrandNode greater than "
	   << 2*nBnd << "." << endl;
      exit(0);
    }
    
    dcn2 = new int*[nStrandNode+1];
    for (int j=0; j<nStrandNode+1; j++) dcn2[j] = new int[2];
    
    dcn2[0][0] = 0; //root stencils
    dcn2[0][1] = 6;
    
    dcn2[1][0] = 0;
    dcn2[1][1] = 6;
    
    dcn2[2][0] = 0;
    dcn2[2][1] = 6;

    dcn2[3][0] = 1;
    dcn2[3][1] = 5;

    dcn2[4][0] = 2;
    dcn2[4][1] = 4;
    
    for (int j=nBnd; j<nStrandNode+1-nBnd; j++){ //interior stencils
      dcn2[j][0] = j-2;
      dcn2[j][1] = 4;
    }
    
    int i=nBnd-1,j1,nj,l;
    for (int j=nStrandNode+1-nBnd; j<nStrandNode+1; j++){ //tip stencils
      j1         = dcn2[i  ][0];
      nj         = dcn2[i--][1];
      dcn2[j][0] = nStrandNode-(j1+nj);
      dcn2[j][1] = nj;
    }
    
    // coefficients
    dcn1 = new double*[nStrandNode+1];
    for (int j=0; j<nStrandNode+1; j++) dcn1[j] = new double[dcn2[j][1]];
    
    dcn1[0][0] = 0.578172459366088; //root coefficients
    dcn1[0][1] =-1.266196663923749;
    dcn1[0][2] = 0.837149798823747;
    dcn1[0][3] =-0.193501484156926;
    dcn1[0][4] = 0.049477930707167;
    dcn1[0][5] =-0.005102040816327;

    dcn1[1][0] =-0.127709893575089;
    dcn1[1][1] = 0.145568041958605;
    dcn1[1][2] = 0.131267445882571;
    dcn1[1][3] =-0.193501484156926;
    dcn1[1][4] = 0.049477930707167;
    dcn1[1][5] =-0.005102040816327;

    dcn1[2][0] = 0.279069767441860;
    dcn1[2][1] =-0.871381110583768;
    dcn1[2][2] = 0.944826767916469;
    dcn1[2][3] =-0.396891314665401;
    dcn1[2][4] = 0.049477930707167;
    dcn1[2][5] =-0.005102040816327;

    dcn1[3][0] = 0.244897959183673;
    dcn1[3][1] =-0.729591836734694;
    dcn1[3][2] = 0.719387755102041;
    dcn1[3][3] =-0.229591836734694;
    dcn1[3][4] =-0.005102040816327;

    dcn1[4][0] = 0.25;
    dcn1[4][1] =-0.75;
    dcn1[4][2] = 0.75;
    dcn1[4][3] =-0.25;
    
    for (int j=nBnd; j<nStrandNode+1-nBnd; j++){ //interior coefficients
      dcn1[j][0] = .25;
      dcn1[j][1] =-.75;
      dcn1[j][2] = .75;
      dcn1[j][3] =-.25;
    }
    
    i = nBnd-1;
    for (int j=nStrandNode+1-nBnd; j<nStrandNode+1; j++){ //tip coefficients
      j1         = dcn2[i][0];
      nj         = dcn2[i][1];
      l          = nj-1;
      for (int m=0; m<nj; m++) dcn1[j][m] =-dcn1[i][l--];
      i--;
    }
  }


  else if (strandOrder == 3){ // Dissipation_6_3
    // dimensions
    int nBnd=7,nMax=9;
    if (nStrandNode < 2*nBnd){
      cout << "\nPlease choose nStrandNode greater than "
	   << 2*nBnd << "." << endl;
      exit(0);
    }
    
    dcn2 = new int*[nStrandNode+1];
    for (int j=0; j<nStrandNode+1; j++) dcn2[j] = new int[2];
    
    dcn2[0][0] = 0; //root stencils
    dcn2[0][1] = 9;
    
    dcn2[1][0] = 0;
    dcn2[1][1] = 9;
    
    dcn2[2][0] = 0;
    dcn2[2][1] = 9;

    dcn2[3][0] = 0;
    dcn2[3][1] = 9;

    dcn2[4][0] = 1;
    dcn2[4][1] = 8;

    dcn2[5][0] = 2;
    dcn2[5][1] = 7;

    dcn2[6][0] = 3;
    dcn2[6][1] = 6;
    
    for (int j=nBnd; j<nStrandNode+1-nBnd; j++){ //interior stencils
      dcn2[j][0] = j-3;
      dcn2[j][1] = 6;
    }
    
    int i=nBnd-1,j1,nj,l;
    for (int j=nStrandNode+1-nBnd; j<nStrandNode+1; j++){ //tip stencils
      j1         = dcn2[i  ][0];
      nj         = dcn2[i--][1];
      dcn2[j][0] = nStrandNode-(j1+nj);
      dcn2[j][1] = nj;
    }
    
    // coefficients
    dcn1 = new double*[nStrandNode+1];
    for (int j=0; j<nStrandNode+1; j++) dcn1[j] = new double[dcn2[j][1]];
    
    dcn1[0][0] = 0.311362984543901; //root coefficients
    dcn1[0][1] =-1.105326803736358;
    dcn1[0][2] = 1.540268773273330;
    dcn1[0][3] =-1.135327578102308;
    dcn1[0][4] = 0.555820097133987;
    dcn1[0][5] =-0.217760746642902;
    dcn1[0][6] = 0.061305150282205;
    dcn1[0][7] =-0.011199448496792;
    dcn1[0][8] = 0.000857571744937;

    dcn1[1][0] = 0.113546294676512;
    dcn1[1][1] =-0.511876734134190;
    dcn1[1][2] = 0.946818703671162;
    dcn1[1][3] =-0.937510888234919;
    dcn1[1][4] = 0.555820097133987;
    dcn1[1][5] =-0.217760746642902;
    dcn1[1][6] = 0.061305150282205;
    dcn1[1][7] =-0.011199448496792;
    dcn1[1][8] = 0.000857571744937;

    dcn1[2][0] = 0.248400202942557;
    dcn1[2][1] =-0.961389761687673;
    dcn1[2][2] = 1.486234336735342;
    dcn1[2][3] =-1.207218704767009;
    dcn1[2][4] = 0.600771399889335;
    dcn1[2][5] =-0.217760746642902;
    dcn1[2][6] = 0.061305150282205;
    dcn1[2][7] =-0.011199448496792;
    dcn1[2][8] = 0.000857571744937;

    dcn1[3][0] =-0.050382534054861;
    dcn1[3][1] = 0.233741186301998;
    dcn1[3][2] =-0.406056330914971;
    dcn1[3][3] = 0.286694980220080;
    dcn1[3][4] = 0.003205925894499;
    dcn1[3][5] =-0.118166500977096;
    dcn1[3][6] = 0.061305150282205;
    dcn1[3][7] =-0.011199448496792;
    dcn1[3][8] = 0.000857571744937;

    dcn1[4][0] =-0.068554018027168;
    dcn1[4][1] = 0.349681679907944;
    dcn1[4][2] =-0.720955700877139;
    dcn1[4][3] = 0.758943936717414;
    dcn1[4][4] =-0.420461705306262;
    dcn1[4][5] = 0.111687684337066;
    dcn1[4][6] =-0.011199448496792;
    dcn1[4][7] = 0.000857571744937;

    dcn1[5][0] =-0.061642428255063;
    dcn1[5][1] = 0.307354569530376;
    dcn1[5][2] =-0.612136423825940;
    dcn1[5][3] = 0.607848565101253;
    dcn1[5][4] =-0.299636423825940;
    dcn1[5][5] = 0.057354569530376;
    dcn1[5][6] = 0.000857571744937;

    dcn1[6][0] =-0.062500000000000;
    dcn1[6][1] = 0.312500000000000;
    dcn1[6][2] =-0.625000000000000;
    dcn1[6][3] = 0.625000000000000;
    dcn1[6][4] =-0.312500000000000;
    dcn1[6][5] = 0.062500000000000;
    
    for (int j=nBnd; j<nStrandNode+1-nBnd; j++){ //interior coefficients
      dcn1[j][0] =-.0625;
      dcn1[j][1] = .3125;
      dcn1[j][2] =-.625;
      dcn1[j][3] = .625;
      dcn1[j][4] =-.3125;
      dcn1[j][5] = .0625;
    }
    
    i = nBnd-1;
    for (int j=nStrandNode+1-nBnd; j<nStrandNode+1; j++){ //tip coefficients
      j1         = dcn2[i][0];
      nj         = dcn2[i][1];
      l          = nj-1;
      for (int m=0; m<nj; m++) dcn1[j][m] =-dcn1[i][l--];
      i--;
    }
  }


  else if (strandOrder == 4){ // Dissipation_8_4
    // dimensions
    int nBnd=9,nMax=12;
    if (nStrandNode < 2*nBnd){
      cout << "\nPlease choose nStrandNode greater than "
	   << 2*nBnd << "." << endl;
      exit(0);
    }
    
    dcn2 = new int*[nStrandNode+1];
    for (int j=0; j<nStrandNode+1; j++) dcn2[j] = new int[2];
    
    dcn2[0][0] = 0; //root stencils
    dcn2[0][1] = 12;
    
    dcn2[1][0] = 0;
    dcn2[1][1] = 12;
    
    dcn2[2][0] = 0;
    dcn2[2][1] = 12;

    dcn2[3][0] = 0;
    dcn2[3][1] = 12;

    dcn2[4][0] = 0;
    dcn2[4][1] = 12;

    dcn2[5][0] = 1;
    dcn2[5][1] = 11;

    dcn2[6][0] = 2;
    dcn2[6][1] = 10;

    dcn2[7][0] = 3;
    dcn2[7][1] = 9;

    dcn2[8][0] = 4;
    dcn2[8][1] = 8;
    
    for (int j=nBnd; j<nStrandNode+1-nBnd; j++){ //interior stencils
      dcn2[j][0] = j-4;
      dcn2[j][1] = 8;
    }
    
    int i=nBnd-1,j1,nj,l;
    for (int j=nStrandNode+1-nBnd; j<nStrandNode+1; j++){ //tip stencils
      j1         = dcn2[i  ][0];
      nj         = dcn2[i--][1];
      dcn2[j][0] = nStrandNode-(j1+nj);
      dcn2[j][1] = nj;
    }
    
    // coefficients
    dcn1 = new double*[nStrandNode+1];
    for (int j=0; j<nStrandNode+1; j++) dcn1[j] = new double[dcn2[j][1]];
    
    dcn1[0][0 ] = 0.379266842169091; //root coefficients
    dcn1[0][1 ] =-1.836668476375949;
    dcn1[0][2 ] = 3.775132781413755;
    dcn1[0][3 ] =-4.440815102468612;
    dcn1[0][4 ] = 3.515177758249977;
    dcn1[0][5 ] =-2.120224880663677;
    dcn1[0][6 ] = 1.012672621119659;
    dcn1[0][7 ] =-0.362917553220777;
    dcn1[0][8 ] = 0.093900382911580;
    dcn1[0][9 ] =-0.017833898057374;
    dcn1[0][10] = 0.002454017821634;
    dcn1[0][11] =-0.000144492899307;

    dcn1[1][0 ] = 0.326281104530594;
    dcn1[1][1 ] =-1.624725525821962;
    dcn1[1][2 ] = 3.457218355582774;
    dcn1[1][3 ] =-4.228872151914625;
    dcn1[1][4 ] = 3.462192020611480;
    dcn1[1][5 ] =-2.120224880663677;
    dcn1[1][6 ] = 1.012672621119659;
    dcn1[1][7 ] =-0.362917553220777;
    dcn1[1][8 ] = 0.093900382911580;
    dcn1[1][9 ] =-0.017833898057374;
    dcn1[1][10] = 0.002454017821634;
    dcn1[1][11] =-0.000144492899307;

    dcn1[2][0 ] = 0.367245353830924;
    dcn1[2][1 ] =-1.798823585348362;
    dcn1[2][2 ] = 3.743968100685081;
    dcn1[2][3 ] =-4.454175523066437;
    dcn1[2][4 ] = 3.544120519212139;
    dcn1[2][5 ] =-2.130465942988759;
    dcn1[2][6 ] = 1.012672621119659;
    dcn1[2][7 ] =-0.362917553220777;
    dcn1[2][8 ] = 0.093900382911580;
    dcn1[2][9 ] =-0.017833898057374;
    dcn1[2][10] = 0.002454017821634;
    dcn1[2][11] =-0.000144492899307;

    dcn1[3][0 ] = 0.003101044789933;
    dcn1[3][1 ] =-0.099483476490407;
    dcn1[3][2 ] = 0.527360037489666;
    dcn1[3][3 ] =-1.298258178044521;
    dcn1[3][4 ] = 1.844780410354184;
    dcn1[3][5 ] =-1.644940197600772;
    dcn1[3][6 ] = 0.951981902946160;
    dcn1[3][7 ] =-0.362917553220777;
    dcn1[3][8 ] = 0.093900382911580;
    dcn1[3][9 ] =-0.017833898057374;
    dcn1[3][10] = 0.002454017821634;
    dcn1[3][11] =-0.000144492899307;

    dcn1[4][0 ] = 0.037859692114567;
    dcn1[4][1 ] =-0.290656036775894;
    dcn1[4][2 ] = 0.979222452709908;
    dcn1[4][3 ] =-1.897844844394457;
    dcn1[4][4 ] = 2.331401472899060;
    dcn1[4][5 ] =-1.888250728873210;
    dcn1[4][6 ] = 1.021499197595428;
    dcn1[4][7 ] =-0.371607215051935;
    dcn1[4][8 ] = 0.093900382911580;
    dcn1[4][9 ] =-0.017833898057374;
    dcn1[4][10] = 0.002454017821634;
    dcn1[4][11] =-0.000144492899307;

    dcn1[5][0 ] = 0.012221500140644;
    dcn1[5][1 ] =-0.080848926497976;
    dcn1[5][2 ] = 0.222297914021312;
    dcn1[5][3 ] =-0.318776975120651;
    dcn1[5][4 ] = 0.231892029542559;
    dcn1[5][5 ] =-0.038572181612456;
    dcn1[5][6 ] =-0.068729678135397;
    dcn1[5][7 ] = 0.056040690797013;
    dcn1[5][8 ] =-0.017833898057374;
    dcn1[5][9 ] = 0.002454017821634;
    dcn1[5][10] =-0.000144492899307;

    dcn1[6][0 ] = 0.016923074627178;
    dcn1[6][1 ] =-0.119904089916727;
    dcn1[6][2 ] = 0.365627032755427;
    dcn1[6][3 ] =-0.623612980302539;
    dcn1[6][4 ] = 0.645831826263622;
    dcn1[6][5 ] =-0.410931682073436;
    dcn1[6][6 ] = 0.153812691922167;
    dcn1[6][7 ] =-0.030055398198018;
    dcn1[6][8 ] = 0.002454017821634;
    dcn1[6][9 ] =-0.000144492899307;

    dcn1[7][0 ] = 0.015480507100693;
    dcn1[7][1 ] =-0.108219056805543;
    dcn1[7][2 ] = 0.324079198819402;
    dcn1[7][3 ] =-0.538783397638803;
    dcn1[7][4 ] = 0.536760497048504;
    dcn1[7][5 ] =-0.320033397638803;
    dcn1[7][6 ] = 0.105329198819402;
    dcn1[7][7 ] =-0.014469056805543;
    dcn1[7][8 ] =-0.000144492899307;

    dcn1[8][0 ] = 0.015625000000000;
    dcn1[8][1 ] =-0.109375000000000;
    dcn1[8][2 ] = 0.328125000000000;
    dcn1[8][3 ] =-0.546875000000000;
    dcn1[8][4 ] = 0.546875000000000;
    dcn1[8][5 ] =-0.328125000000000;
    dcn1[8][6 ] = 0.109375000000000;
    dcn1[8][7 ] =-0.015625000000000;
    
    for (int j=nBnd; j<nStrandNode+1-nBnd; j++){ //interior coefficients
      dcn1[j][0] = 0.015625;
      dcn1[j][1] =-0.109375;
      dcn1[j][2] = 0.328125;
      dcn1[j][3] =-0.546875;
      dcn1[j][4] = 0.546875;
      dcn1[j][5] =-0.328125;
      dcn1[j][6] = 0.109375;
      dcn1[j][7] =-0.015625;
    }
    
    i = nBnd-1;
    for (int j=nStrandNode+1-nBnd; j<nStrandNode+1; j++){ //tip coefficients
      j1         = dcn2[i][0];
      nj         = dcn2[i][1];
      l          = nj-1;
      for (int m=0; m<nj; m++) dcn1[j][m] =-dcn1[i][l--];
      i--;
    }
  }


  else{
    cout << "\nThis strand order not yet available for dissipation." << endl;
    exit(0);
  }

  /*
  int nj,j1;
  for (int j=0; j<nStrandNode+1; j++){
    j1 = dcn2[j][0];
    nj = dcn2[j][1];
    cout << "\nEdge: " << j << " " << endl;
    for (int m=0; m<nj; m++) cout << j1++ << " " << dcn1[j][m] << endl;
  }
  exit(0);
  */
}
