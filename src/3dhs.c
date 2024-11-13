/**
 * Author: Aman Hogan-Bailey
 * Solves 3d heat equation for one mpi proccess
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include "../include/3dsolver.h"

int main(int argc, char **argv) 
{
    // Initialize MPI
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Check we are using one mpi proccess
    if (size != 1) 
    {
        if (rank == 0) 
        {
            fprintf(stderr, "Error: This program should be run with exactly one MPI process.\n");
        }
        MPI_Finalize();
        return -1;
    }

    // Final time
    double T = 1.0; 

    // Grid sizes for 3D
    int GRID_SIZES[] = {25,50,100}; 
    
    // CFL condition
    double CFL = 0.5; 

    // For each grid size, calcualte error and runtime
    for (int k = 0; k < 3; k++) 
    {
        // Initialize heat equation variables
        int GRID_SIZE = GRID_SIZES[k];
        double dx = PI / GRID_SIZE;
        double dy = PI / GRID_SIZE;
        double dz = PI / GRID_SIZE;
        double dx2 = dx * dx;
        double dy2 = dy * dy;
        double dz2 = dz * dz;
        double dt = CFL * dx2 / 3.0;
        int nsteps = (int)ceil(T / dt);
        dt = T / nsteps;

        double *u = (double *)malloc((GRID_SIZE + 1) * (GRID_SIZE + 1) * (GRID_SIZE + 1) * sizeof(double));
        double *u_new = (double *)malloc((GRID_SIZE + 1) * (GRID_SIZE + 1) * (GRID_SIZE + 1) * sizeof(double));
        
        initialize(u, GRID_SIZE, dx, dy, dz);
        double start_time = MPI_Wtime();

        for (int t_step = 0; t_step < nsteps; t_step++)
        {
            time_step(u, u_new, GRID_SIZE, dt, dx2, dy2, dz2);
            double *tmp = u;
            u = u_new;
            u_new = tmp;
        }
        
        double end_time = MPI_Wtime();
        double error = compute_error(u, GRID_SIZE, T, dx, dy, dz);
        double runtime = end_time - start_time;
        printf("Grid size N = %d, Error = %lf, Runtime = %lf seconds\n", GRID_SIZE, error, runtime);

        free(u);
        free(u_new);
    }

    MPI_Finalize();
    return 0;
}
