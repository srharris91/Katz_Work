#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define UP      0
#define DOWN    1
#define LEFT    2
#define RIGHT   3

using namespace std;

int main (int argc, char** argv)
	{
	int N = 120;
    int rank, numtasks, numtasks_all, new_numtasks, rank_old;
    int reorder = 1;
    int coords[2];
    int nbrs[4];
    MPI_Comm GRID_COMM, new_comm;
    
    int pbc_check[2]    = {0,0};

    MPI_Init(&argc, &argv);
   

    MPI_Comm_size(MPI_COMM_WORLD, &numtasks_all);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_old);

    int size_grid = int(floor(sqrt(double(numtasks_all))))*int(floor(sqrt(double(numtasks_all))));
    int size_grid_ranks[size_grid];

    //cout << size_grid << endl << endl;

    for(int i = 0; i<size_grid; i++){size_grid_ranks[i] = i;}

    MPI_Group old_group, new_group;
    MPI_Comm_group(MPI_COMM_WORLD, &old_group);

    MPI_Group_incl(old_group, size_grid, size_grid_ranks, &new_group);
    MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm);

    if (rank_old >= size_grid){
        cout << "**********Rank*********** = " << rank_old << endl;
        MPI_Finalize();
        return 0;}

    MPI_Comm_size(new_comm, &new_numtasks);

    cout << "Old Rank: " << rank_old << "\tNew Rank: " << rank << endl;
    
    //int dims[2] = {size_grid,size_grid};
    int dims[2] = {0,0};
    MPI_Barrier(new_comm);
    MPI_Dims_create(new_numtasks, 2, dims);

    MPI_Cart_create(new_comm, 2, dims, pbc_check, reorder, &GRID_COMM);

    MPI_Comm_free(&new_comm);

    MPI_Comm_size(GRID_COMM, &numtasks);
    MPI_Comm_rank(GRID_COMM, &rank);


    
    if(rank == 0){cout << "dims = (" << dims[1] << "," << dims[0] << ")" << endl << endl;}
    
    int size_h = N/dims[1] + 2;
    int size_v = N/dims[0] + 2;

    MPI_Cart_coords(GRID_COMM, rank, 2, coords);
    MPI_Cart_shift(GRID_COMM, 0, 1, &nbrs[UP],    &nbrs[DOWN]);
    MPI_Cart_shift(GRID_COMM, 1, 1, &nbrs[LEFT],  &nbrs[RIGHT]);

    //cout <<"rank : " << rank << "\t Coordinates : (" << coords[0] << ", " << coords[1] << ")" << endl;

    printf("My rank is %d      My neighbors are (UP = %d, DOWN = %d, LEFT = %d, RIGHT = %d)\n", rank, nbrs[UP], nbrs[DOWN], nbrs[LEFT], nbrs[RIGHT]);
    //printf("My rank is %d       My coordinates are(%d, %d) \n", rank, coords[0], coords[1]);


    // Creation of MPI Struct that will be passed between processors
    MPI_Datatype mpi_p;
    int count = 1;
    int blens[1]        = {2};
    MPI_Aint index[1]   = {0};
    MPI_Datatype old_types[1] = {MPI_DOUBLE};

    MPI_Type_struct(count, blens, index, old_types, &mpi_p);
    MPI_Type_commit(&mpi_p);


    //MPI_Request req_v[4*size_h-8],   req_h[4*size_v-8];
    //MPI_Status  stats_v[4*size_h-8], stats_h[4*size_v-8];

    struct rs{double p_1, p_old;};

	rs phi[size_h][size_v];
    //double phi_old  [N/dims[0]+2][N/dims[1]+2];

		
	double Beta = 0.99;
	//cout << "What is the Blending Factor? (0 <= Beta <= 1)" << endl;
	//cin >> Beta;
	//cout << endl;
	
//	cout << "How many iterations would you like to run?" << endl;
//	cin >> it;
    int it = 1;
	
	// Set phi equal to 0
	for (int i = 0; i<size_h; i++){
		for ( int j = 0; j<size_v; j++ )
			{phi[i][j].p_1   = 0.;
			 phi[i][j].p_old = 0.;}}
	
	// Set Boundary Conditions on the West and South Faces
	for (int j = 1; j<size_v-1; j++){
		// West Side
            if (coords[1] == 0) { 
			    phi[0][j].p_1   = 100.;
			    phi[0][j].p_old = 100.;}}
    
    //if(rank == 0){
        //for(int j = 1; j<size_v-1; j++){cout << "phi[0][" << j << "] = " << phi[0][j].p_1 << endl;}}

    for (int i = 1; i<size_h-1; i++){
		// South Side
            if (coords[0] == dims[0]-1) {
			    phi[i][size_v+1].p_1   = 0.;
			    phi[i][size_v+1].p_old = 0.;}}
    
    //if (rank == 0){cout << "Dims = (" << dims[0] << ", " << dims[1] << ")" << endl << endl;}

//Upwind Coefficients
double a_W = 2., a_S = 2., a_P = a_W + a_S;

//Central Difference Coefficients
double a_cE = -1., a_cN = -1., a_cW = 1., a_cS = 1., a_cP = a_cE + a_cN + a_cW + a_cS;

//Initialize error variables
double e, res = 1., new_res = 1.;

// Residuals File	
ofstream residual;
residual.open ("Figures/res/resid_40_mpi.txt");

//Computation of phi
//for ( int n = 1; n<it; n++ )

int n = 1, r = 0;
double tol = 1e-10;
while (new_res > tol){
	e = 0., res = 0., new_res = 0.;
	for (int i = 1; i<size_h-2; i++){
			for (int j = 1; j<size_v-2; j++){
					if (coords[1] == 0 && i == 1)
					{ // West Boundary Condition
					phi[i][j].p_1 = (a_W*phi[i-1][j].p_1 + a_S*phi[i][j-1].p_1)/a_P - (Beta/a_P)*(a_cP*phi[i][j].p_old -
						a_cE*phi[i+1][j].p_old - 200. - a_cN*phi[i][j+1].p_old - a_cS*phi[i][j-1].p_old - 
						a_P*phi[i][j].p_old + a_W*phi[i-1][j].p_old + a_S*phi[i][j-1].p_old);
					}
					if (coords[1] == dims[1]-1  && i == size_h-3)  // may have to change to dims[0]-1
					{ // East Boundary Condition
					phi[i][j].p_1 = (a_W*phi[i-1][j].p_1 + a_S*phi[i][j-1].p_1)/a_P - (Beta/a_P)*(a_cP*phi[i][j].p_old -
						a_cE*phi[i][j].p_old - a_cW*phi[i-1][j].p_old - a_cN*phi[i][j+1].p_old - a_cS*phi[i][j-1].p_old - 
						a_P*phi[i][j].p_old + a_W*phi[i-1][j].p_old + a_S*phi[i][j-1].p_old);
					}
					else if (coords[0] == 0 && j == 1)
					{ // Top Boundary Condition
					phi[i][j].p_1 = (a_W*phi[i-1][j].p_1 + a_S*phi[i][j-1].p_1)/a_P - (Beta/a_P)*(a_cP*phi[i][j].p_old -
						a_cE*phi[i+1][j].p_old - a_cW*phi[i-1][j].p_old - a_cN*phi[i][j].p_old - a_cS*phi[i][j-1].p_old - 
						a_P*phi[i][j].p_old + a_W*phi[i-1][j].p_old + a_S*phi[i][j-1].p_old);
					}
					else if (coords[0] == dims[0]-1 && j == size_v-3)
					{ // Bottom Boundary Condition
					phi[i][j].p_1 = (a_W*phi[i-1][j].p_1 + a_S*phi[i][j-1].p_1)/a_P - (Beta/a_P)*(a_cP*phi[i][j].p_old -
						a_cE*phi[i+1][j].p_old - a_cW*phi[i-1][j].p_old - a_cN*phi[i][j+1].p_old - 0. - 
						a_P*phi[i][j].p_old + a_W*phi[i-1][j].p_old + a_S*phi[i][j-1].p_old);
					}
					else
					{
					phi[i][j].p_1 = (a_W*phi[i-1][j].p_1 + a_S*phi[i][j-1].p_1)/a_P - (Beta/a_P)*(a_cP*phi[i][j].p_old -
						a_cE*phi[i+1][j].p_old - a_cW*phi[i-1][j].p_old - a_cN*phi[i][j+1].p_old - a_cS*phi[i][j-1].p_old - 
						a_P*phi[i][j].p_old + a_W*phi[i-1][j].p_old + a_S*phi[i][j-1].p_old);
					}
                    
                    if(rank == 2 && i == 1 && (it > 10 && it < 15)){
                        //cout << "phi[" << i << "][" << j << "] = " << phi[i][j].p_1 << endl;
                        }

					e = phi[i][j].p_old - phi[i][j].p_1;
			
					res += e*e;	
				}

		}

        
        //MPI_Barrier(GRID_COMM);
        //if (rank == 0 && (it > 10 && it < 15)){
        //cout << "===============================" << endl;}

        //if (it == 1){cout << "Rank = " << rank << "\tCoordinates = (" << coords[0] << ", " << coords[1] << ")" << endl;
                     //cout << "My Neighbors are: UP = " << nbrs[UP] << "  DOWN = " << nbrs[DOWN] << "  RIGHT = " << nbrs[RIGHT] << "  LEFT = " << nbrs[LEFT] << endl << endl;}


        //MPI_Barrier(GRID_COMM);
        //cout << "P" << rank << endl;
        //MPI_Barrier(GRID_COMM);

        MPI_Allreduce(&res, &new_res, 1, MPI_DOUBLE, MPI_SUM, GRID_COMM);
        //MPI_Barrier(GRID_COMM);

        new_res = sqrt(new_res/(N*N));

    if(it > 5e3){break;}

		for ( int i = 0; i<size_h; i++ )
		{for ( int j = 0; j<size_v; j++ ){phi[i][j].p_old = phi[i][j].p_1;}}
        

        MPI_Barrier(GRID_COMM);
        if(rank == 0 
        && it % 3 == 0
        ){
        //cout << "Iteration: " << it << "\t Residual: " << new_res << endl;

        residual << it << "\t" << new_res << endl;
        }

        it++;
	    
        MPI_Request req_v[4*(size_h-2)],   req_h[4*(size_v-2)];
        MPI_Status  stats_v[4*(size_h-2)], stats_h[4*(size_v-2)];

    if(dims[0] > 1){ 
        for (int i = 1; i<size_h-1; i++){
            MPI_Isend(&phi[i][1],        1, mpi_p, nbrs[UP],    0, GRID_COMM, &req_v[i-1]);                 // Send to up
            MPI_Isend(&phi[i][size_h-2], 1, mpi_p, nbrs[DOWN],  0, GRID_COMM, &req_v[i-1+(size_h-2)]);      // Send to down

            MPI_Irecv(&phi[i][0],        1, mpi_p, nbrs[UP],    0, GRID_COMM, &req_v[i-1+2*(size_h-2)]);    // Receive from up
            MPI_Irecv(&phi[i][size_h-1], 1, mpi_p, nbrs[DOWN],  0, GRID_COMM, &req_v[i-1+3*(size_h-2)]);    // Receive from down
        }
        
        MPI_Waitall(4*(size_h-2), req_v, stats_v);}

    if(dims[1] > 1){
        for(int j = 1; j<size_v-1; j++){
            MPI_Isend(&phi[1][j],        1, mpi_p, nbrs[LEFT],  0, GRID_COMM, &req_h[j-1]);                 // Send to left
            MPI_Isend(&phi[size_v-2][j], 1, mpi_p, nbrs[RIGHT], 0, GRID_COMM, &req_h[j-1+(size_v-2)]);      // Send to right

            MPI_Irecv(&phi[0][j],        1, mpi_p, nbrs[LEFT],  0, GRID_COMM, &req_h[j-1+2*(size_v-2)]);    // Receive from left
            MPI_Irecv(&phi[size_v-1][j], 1, mpi_p, nbrs[RIGHT], 0, GRID_COMM, &req_h[j-1+3*(size_v-2)]);    // Receive from right
        }

        MPI_Waitall(4*(size_v-2), req_h, stats_h);}
        
        //for(int i = 0; i<4*size_h-8; i++){MPI_Request_free(&req_v[i]);}
        //MPI_Grequest_complete(req_v);

        //for(int j = 0; j<4*size_v-8; j++){MPI_Request_free(&req_h[j]);}
        //MPI_Grequest_complete(req_h);
        
		if (r == 100){
            r = 0;}
		//	cout << n << "\t" << new_res << endl;}

		//for ( int i = 0; i<size_h; i++ )
		//{for ( int j = 0; j<size_v; j++ ){phi[i][j].p_old = phi[i][j].p_1;}}
		n++;
		r++;
	}



residual.close();


MPI_Finalize();

if(rank == 0){
	cout << "Number of Iterations was: " << n-1 << "\tFinal Residual: " << new_res << endl;}
/*
// Contour Plot file
	ofstream file;
	file.open ("results_40_1.txt");

	for ( int i = 0; i<N+2; i++ )
	{
		for ( int j = 0; j<N+2; j++ )
		{
			if ( i == 0 )
			{
				if ( j == 0 )
				{
					
				file << 0 << "\t" 
					<< 0 << "\t" 
					<< phi[i][j] << endl;
				}
				else if ( j == N+1 )
				{
				file << 0 << "\t" << 1 << "\t" << phi[i][j-1] << endl;
				}
				else
				file << 0 << "\t" 
					 << double(j)/double(N) - 1./(2.*double(N)) << "\t"
					 << phi[i][j] << endl;
			}
			
			else if ( i == N+1 )
			{
			
			if ( j == 0 )
				file << 1 << "\t" << 0 << "\t" << 0 << endl;
			else if ( j == N+1 )
				file << 1 << "\t" << 1 << "\t" << phi[i][j] << endl;
			
			else
				file << 1 << "\t" 
					<< double(j)/double(N) - 1./(2.*double(N)) << "\t" 
					<< phi[i-1][j] << endl;
			}
			
			else if ( j == 0 )
			{
			file << double(i)/double(N) - 1./(2.*double(N)) << "\t"<< 0 << "\t" 
				 << phi[i][j] << endl;
			}
			
			else if ( j == N+1 )
			{
			file << double(i)/double(N)  - 1./(2.*double(N)) << "\t" << 1 << "\t" 
				 << phi[i][j-1]<< endl;
			}
			else 
			{
			file << double(i)/double(N) - 1./(2.*double(N)) << "\t" 
				 << double(j)/double(N) - 1./(2.*double(N)) << "\t" 
				 << phi[i][j] << endl;
			}
				 
		};
		file << endl;
	};
	
	file.close();
	
// Error Plot File
	ofstream error;
	error.open ("error_40_.txt");
	

	for ( int i = 0; i<N+2; i++ )
	{
		if ( i == 0 ) 
			error << 0 << "\t" << phi[i][(N+1)-i] << endl;
		else if ( i == N+1 )
			error << 1 << "\t" << phi[i][(N+1)-1] << endl;
		else
			error << double(i)/double(N) - 1./(2.*double(N)) << "\t" << phi[i][(N+1)-i] << endl;
	};
	
	error.close();
*/	
//if(rank == 0){system("pdflatex Figures/res/residuals.tex");}

	return 0;
	}
