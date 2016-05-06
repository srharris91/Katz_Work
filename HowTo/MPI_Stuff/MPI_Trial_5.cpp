/*
<!-- saved from url=(0049)http://grid.hpctools.uh.edu/6397/slides/MPI/cfd.c -->
*******************************************************************************
The author of this code is Noah Yan yanyh@cs.uh.edu
And it is used for the course testing and it is not suitable for any production usage. If you have any questions, free free to ask the author
*******************************************************************************
*/

#include "mpi.h"
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

// #define num_grid 10000       // Number of Grid configured to compute
#define itmax 20000     // max loop times
#define stride 1            // After "strid" loop, to syn to calculate the termination condition
#define true    1
#define false   0
#define max(a, b) ((a) > (b) ? (a) : (b))
#define addr(xx, yy, zz) ((xx)*(ty+2)*(tz+2) + (yy)*(tz+2) + (zz))

// These are used to store the rank of the neighbour
int right, left, front, back, upper, lower;
int ierror, my_rank, size;
int cord[3];

const double edvu = 1.0e-05;
MPI_Comm cube;

int main(int argc, char *argv[])
{
    int tx_from, ty_from, tz_from, tx, ty, tz, qt_x, qt_y, qt_z, num_grid;
    double *tp, *tpt, *tptr, tmp, tmptp;
    double dx, dt, dx2, diff, maxdf, smaxdf;
    double std;
    double starttime, endtime;
    double commtime;
    int i, j, k, it, num_sub, cal_grid;
    MPI_Datatype xvtr, yvtr, zvtr;
    MPI_Status status[4];
    MPI_Request rqst[4];
    int dim[3], periods[3];
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    num_grid=atoi(argv[1]);
    std = 1.0 / num_grid;
    // Get the information about how to split the cube
    dim[0]=dim[1]=dim[2]=0;
    MPI_Dims_create(size, 3, dim);

    // Create the virtual topology which is a 3D cube. Each part of the
    // cube donates those points which should be computed by my_rank processes
    // my_rank process is also mapped to its cordinate in three dimentions represented by cord[3] array
    periods[0] = false;
    periods[1] = false;
    periods[2] = false;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dim, periods, false, &cube);
    MPI_Cart_coords(cube, my_rank, 3, cord);    // cord[0]=x, cord[1]=y, cord[2]=z

    // Calcuate the bounds of those points which are computed by my_rank(now is mapped to cord[])
    //tx_from is the starting calculation point in x-cord, tx is the number of points for this processes. So is ty, tz.
    cal_grid = num_grid - 1;
    qt_x = cal_grid / dim[0];
    front = dim[0] - cal_grid % dim[0];
    tx_from = 1 + cord[0] * qt_x + max(0, cord[0] - front);
    tx = qt_x;
    if ((cord[0] - front) >= 0) tx++;

    qt_y = cal_grid / dim[1];
    front = dim[1] - cal_grid % dim[1];
    ty_from = 1 + cord[1] * qt_y + max(0, cord[1] - front);
    ty = qt_y;
    if ((cord[1] - front) >= 0) ty++;

    qt_z = cal_grid / dim[2];
    front = dim[2] - cal_grid % dim[2];
    tz_from = 1 + cord[2] * qt_z + max(0, cord[2] - front);
    tz = qt_z;
    if ((cord[2] - front) >= 0) tz++;

    // Create the dynamic storage for computing
    num_sub = (tx + 2) * (ty + 2) * (tz + 2);
    tp = (double *) malloc(sizeof(double) * num_sub);
    tpt = (double *) malloc(sizeof(double) * num_sub);

    // Itialization of the array to 0, including the boundary data
    for (i = 0; i < num_sub; i++) {
        tp[i] = 0.0;
        tpt[i] = 0.0;
    }

    // Find the neighbour of my_rank by calling MPI_Cart_shift(), and
    // stores the neighbours in the right, left, front, back, upper and lower.MPI_PROC_NULL is 0?
    // Together to itialize those neighbour data for sending and receiving

    MPI_Cart_shift(cube, 0, 1, &left, &right);
    MPI_Cart_shift(cube, 1, 1, &back, &front);
    MPI_Cart_shift(cube, 2, 1, &upper, &lower);

    if (right == MPI_PROC_NULL)
        for (j = 0; j < ty + 2; j++)
            for (k = 0; k < tz + 2; k++)
                tp[addr(tx + 1, j, k)] = 1;

    if (front == MPI_PROC_NULL)
        for (i = 0; i < tx + 2; i++)
            for (k = 0; k < tz + 2; k++)
                tp[addr(i, ty + 1, k)] = (tx_from + i - 1) * std;

    if (back == MPI_PROC_NULL)
        for (i = 0; i < tx + 2; i++)
            for (k = 0; k < tz + 2; k++)
                tp[addr(i, 0, k)] = (tx_from + i - 1) * std;

    if (upper == MPI_PROC_NULL)
        for (i = 0; i < tx + 2; i++)
            for (j = 0; j < ty + 2; j++)
                tp[addr(i, j, tz + 1)] = (tx_from + i - 1) * std;

    if (lower == MPI_PROC_NULL)
        for (i = 0; i < tx + 2; i++)
            for (j = 0; j < ty + 2; j++)
                tp[addr(i, j, 0)] = (tx_from + i - 1) * std;

    // End of Intialization
    // Create the MPI derived  vector type for data exchange in three dimentional
    MPI_Type_contiguous((ty + 2) * (tz + 2), MPI_DOUBLE, &xvtr);
    MPI_Type_commit(&xvtr);
    MPI_Type_vector(tx + 2, tz + 2, (ty + 2) * (tz + 2), MPI_DOUBLE, &yvtr);
    MPI_Type_commit(&yvtr);
    MPI_Type_vector((tx + 2) * (ty + 2), 1, tz + 2, MPI_DOUBLE, &zvtr);
    MPI_Type_commit(&zvtr);

    // Calculating, dx=dy=dz=std, dt=(1/8)*dx2
    dx = std;
    dx2 = dx * dx;
    dt = 1.0 / 8.0 * dx2;
    MPI_Barrier(cube);      // Syn all processes

    starttime = MPI_Wtime();
    commtime = 0;

    for (it = 1; it < itmax; it++) {    // Main loop
        maxdf = 0.0;
        for (i = 1; i < tx + 1; i++)
            for (j = 1; j < ty + 1; j++)
                for (k = 1; k < tz + 1; k++) {
                    tmp = tp[addr(i, j, k)];
                    tmptp =
                                    tp[addr(i + 1, j, k)] + tp[addr(i - 1, j, k)] - 2.0 * tmp +
                                    tp[addr(i, j + 1, k)] + tp[addr(i, j - 1, k)] - 2.0 * tmp +
                                    tp[addr(i, j, k + 1)] + tp[addr(i, j, k - 1)] - 2.0 * tmp;
                    tmptp = tmptp * dt / (2.0 * dx2);
                    tpt[addr(i, j, k)] = tmptp + tmp;
                    diff = tmptp;
                    maxdf = max(diff, maxdf);
                }
    
    // Exchage the old and new data
    for (i = 0; i < num_sub; i++) tp[i] = tpt[i];

    if(( it % stride ) == 0){
        smaxdf = maxdf;
        MPI_Allreduce(&smaxdf, &maxdf, 1, MPI_DOUBLE, MPI_MAX, cube);
        if (maxdf < edvu)                       // end of all the computation
            goto loopend;
    }   
    
    // Send data and receive data to and from neighbour
    commtime = commtime - MPI_Wtime();
    i = 0;
    if (right != MPI_PROC_NULL){
        ierror=MPI_Irecv(&tp[addr(tx + 1, 0, 0)], 1, xvtr, right, MPI_ANY_TAG, cube, &rqst[i++]);
        ierror=MPI_Isend(&tp[addr(tx, 0, 0)], 1, xvtr, right, 0, cube, &rqst[i++]);
    }

    if (left != MPI_PROC_NULL){
        ierror=MPI_Irecv(&tp[addr(0, 0, 0)], 1, xvtr, left, MPI_ANY_TAG, cube, &rqst[i++]);
        ierror=MPI_Isend(&tp[addr(1, 0, 0)], 1, xvtr, left, 0, cube, &rqst[i++]);
    }
    MPI_Waitall(i, rqst, status);

    i = 0;
    if (front != MPI_PROC_NULL){
        ierror=MPI_Irecv(&tp[addr(0, ty + 1, 0)], 1, yvtr, front, MPI_ANY_TAG, cube, &rqst[i++]);
        ierror=MPI_Isend(&tp[addr(0, ty, 0)], 1, yvtr, front, 0, cube, &rqst[i++]);
    }

    if (back != MPI_PROC_NULL){
        ierror=MPI_Irecv(&tp[addr(0, 0, 0)], 1, yvtr, back, MPI_ANY_TAG, cube, &rqst[i++]);
        ierror=MPI_Isend(&tp[addr(0, 1, 0)], 1, yvtr, back, 0, cube, &rqst[i++]);
    }
    MPI_Waitall(i, rqst, status);

    i = 0;
    if (upper != MPI_PROC_NULL){
        ierror=MPI_Irecv(&tp[addr(0, 0, tz + 1)], 1, zvtr, upper, MPI_ANY_TAG, cube, &rqst[i++]);
        ierror=MPI_Isend(&tp[addr(0, 0, tz)], 1, zvtr, upper, 0, cube, &rqst[i++]);
    }

    if (lower != MPI_PROC_NULL){
        ierror=MPI_Irecv(&tp[addr(0, 0, 0)], 1, zvtr, lower, MPI_ANY_TAG, cube, &rqst[i++]);
        ierror=MPI_Isend(&tp[addr(0, 0, 1)], 1, zvtr, lower, 0, cube, &rqst[i++]);
    }
    MPI_Waitall(i, rqst, status);
    
    commtime = commtime + MPI_Wtime();
    }
    loopend:
    it--;
    endtime = MPI_Wtime();
    if (my_rank==0){
        printf("\nNumber of Grids: %d", num_grid);
        printf("\nsize     iter-      Total loop time        comm time     ");
        printf("\n         ations       [seconds]            [seconds]     ");
        printf("\n%4d  %6d        %3.12f         %3.9f       \n",size, it, endtime - starttime, commtime);
    }
    MPI_Finalize();
}

