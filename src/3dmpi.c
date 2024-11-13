/**
 * Author: Aman Hogan-Bailey
 * Solves 3d heat equation using mpi
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include "../include/3dsolver.h"

#define GRID_SIZE 120
#define TIME_FINAL 1.0

int main(int argc, char **argv) 
{
    // Intialize MPI environment
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int MPIX = 2; // Num of chunks X axis of grid is split into
    int MPIY = 2;  // Num of chunks Y axis of grid is split into
    int MPIZ = 2; // Num of chunks Z axis of grid is split into

    if (argc > 1) MPIX = atoi(argv[1]);
    if (argc > 2) MPIY = atoi(argv[2]);
    if (argc > 3) MPIZ = atoi(argv[3]);
    
    int nx = GRID_SIZE / MPIX; // Num of x points per proccess
    int ny = GRID_SIZE / MPIY; // Num of y points per proccess
    int nz = GRID_SIZE / MPIZ; // Num of z points per proccess

    printf("INFO_START: Using %d proccesses\n", MPIX*MPIY*MPIZ);

    // Check if the nprocs match chunks
    if (size != MPIX * MPIY * MPIZ) 
    {
        if (rank == 0) 
        {
            fprintf(stderr, "ERROR: Nprocs (%d) != Grid (%dx%dx%d=%d).\n",size,MPIX,MPIY,MPIZ,MPIX*MPIY*MPIZ);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // Set up commincation for cartiseian topology
    int dims[3] = {MPIX, MPIY, MPIZ};
    int periods[3] = {0, 0, 0};
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &cart_comm);
    int coords[3];
    MPI_Cart_coords(cart_comm, rank, 3, coords);

    // Intialize heat equation variables
    double dx = PI / GRID_SIZE;
    double dy = PI / GRID_SIZE;
    double dz = PI / GRID_SIZE;
    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double dz2 = dz * dz;
    double CFL = 0.5;
    double dt = CFL * dx2 / 3.0;
    int nsteps = (int)ceil(TIME_FINAL / dt);
    dt = TIME_FINAL / nsteps;
    double *u = (double *)malloc((nx + 2) * (ny + 2) * (nz + 2) * sizeof(double));
    double *u_new = (double *)malloc((nx + 2) * (ny + 2) * (nz + 2) * sizeof(double));

    // Intialize ghost cell offsets
    int x_offset = coords[0] * nx;
    int y_offset = coords[1] * ny;
    int z_offset = coords[2] * nz;

    initialize_mpi(u, nx, ny, nz, dx, dy, dz, x_offset, y_offset, z_offset);

    // Vector to store ranks of neighboring proccesses
    int neighbors[6];

    // Assign ranks to proccesses
    MPI_Cart_shift(cart_comm, 2, 1, &neighbors[0], &neighbors[1]);
    MPI_Cart_shift(cart_comm, 1, 1, &neighbors[2], &neighbors[3]);
    MPI_Cart_shift(cart_comm, 0, 1, &neighbors[4], &neighbors[5]);

    
    printf("INFO: Process %d: (x, y, z) = (%d, %d, %d) | Neighbors - Bottom: %d, Top: %d, Left: %d, Right: %d, Front: %d, Back: %d\n",
    rank, coords[0], coords[1], coords[2], neighbors[0], neighbors[1], neighbors[2], neighbors[3], neighbors[4], neighbors[5]);

    double start_time, end_time;
    start_time = MPI_Wtime();

    // For each given timestep, update u
    for (int t = 0; t < nsteps; t++) 
    {
        time_step_mpi(u, u_new, nx, ny, nz, dt, dx2, dy2, dz2, cart_comm, neighbors);
        double *tmp = u;
        u = u_new;
        u_new = tmp;
    }

    end_time = MPI_Wtime();

    // Get local and global error of processes
    double local_error = compute_local_error_mpi(u, nx, ny, nz, TIME_FINAL, dx, dy, dz, x_offset, y_offset, z_offset);
    double global_error_squared;
    MPI_Reduce(&local_error, &global_error_squared, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) 
    {
        double global_error = sqrt(global_error_squared);
        printf("GLOBAL_ERROR: %f\n", global_error);
        printf("RUNTIME: %f seconds\n", end_time - start_time);
    }

    free(u);
    free(u_new);
    MPI_Finalize();
    return 0;
}
