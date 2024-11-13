/**
 * Author: Aman Hogan-Bailey
 * functions file that contains functions for 
 * 3d heatsolver for both mpi and mpi on one proccess
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "../include/3dsolver.h"

void initialize_mpi(double *u, int nx, int ny, int nz, double dx, double dy, double dz, int xo, int yo, int zo) 
{
    for (int i = 1; i <= nx; i++) 
    {
        for (int j = 1; j <= ny; j++) 
        {
            for (int k = 1; k <= nz; k++) 
            {
                double x = (i - 1 + xo) * dx;
                double y = (j - 1 + yo) * dy;
                double z = (k - 1 + zo) * dz;
                u[I3DMPI(i, j, k, nx, ny)] = sin(x) * sin(y) * sin(z);
            }
        }
    }
}


void time_step_mpi(double *u, double *u_new, int nx, int ny, int nz, double dt, double dx2, double dy2, double dz2, MPI_Comm cart_comm, int neighbors[6]) 
{
    int rank;
    MPI_Status status;
    MPI_Comm_rank(cart_comm, &rank);

    // Z TOP & BOTTOM
    if (neighbors[0] != MPI_PROC_NULL || neighbors[1] != MPI_PROC_NULL) 
    {
        MPI_Request requests[4];
        MPI_Status statuses[4];
        int req_count = 0;

        // RECIEVES
        if (neighbors[0] != MPI_PROC_NULL) { MPI_Irecv(&u[I3DMPI(1, 1, 0, nx, ny)], nx * ny, MPI_DOUBLE, neighbors[0], 1, cart_comm, &requests[req_count++]);}
        if (neighbors[1] != MPI_PROC_NULL) { MPI_Irecv(&u[I3DMPI(1, 1, nz + 1, nx, ny)], nx * ny, MPI_DOUBLE, neighbors[1], 0, cart_comm, &requests[req_count++]);}

        // SENDS
        if (neighbors[0] != MPI_PROC_NULL) { MPI_Isend(&u[I3DMPI(1, 1, 1, nx, ny)], nx * ny, MPI_DOUBLE, neighbors[0], 0, cart_comm, &requests[req_count++]); }
        if (neighbors[1] != MPI_PROC_NULL) { MPI_Isend(&u[I3DMPI(1, 1, nz, nx, ny)], nx * ny, MPI_DOUBLE, neighbors[1], 1, cart_comm, &requests[req_count++]); }

        MPI_Waitall(req_count, requests, statuses);
    }

    // Y LEFT & RIGHT
    if (neighbors[2] != MPI_PROC_NULL || neighbors[3] != MPI_PROC_NULL) 
    {
        MPI_Request requests[4];
        MPI_Status statuses[4];
        int req_count = 0;

        // RECIEVES
        if (neighbors[2] != MPI_PROC_NULL) { MPI_Irecv(&u[I3DMPI(1, 0, 1, nx, ny)], nx * nz, MPI_DOUBLE, neighbors[2], 3, cart_comm, &requests[req_count++]);}
        if (neighbors[3] != MPI_PROC_NULL) { MPI_Irecv(&u[I3DMPI(1, ny + 1, 1, nx, ny)], nx * nz, MPI_DOUBLE, neighbors[3], 2, cart_comm, &requests[req_count++]);}

        // SENDS
        if (neighbors[2] != MPI_PROC_NULL) {MPI_Isend(&u[I3DMPI(1, 1, 1, nx, ny)], nx * nz, MPI_DOUBLE, neighbors[2], 2, cart_comm, &requests[req_count++]);}
        if (neighbors[3] != MPI_PROC_NULL) {MPI_Isend(&u[I3DMPI(1, ny, 1, nx, ny)], nx * nz, MPI_DOUBLE, neighbors[3], 3, cart_comm, &requests[req_count++]);}

        MPI_Waitall(req_count, requests, statuses);
    }

    // X FRONT & BACK
    if (neighbors[4] != MPI_PROC_NULL || neighbors[5] != MPI_PROC_NULL) 
    {
        MPI_Request requests[4];
        MPI_Status statuses[4];
        int req_count = 0;

        // RECEIVES
        if (neighbors[4] != MPI_PROC_NULL) { MPI_Irecv(&u[I3DMPI(0, 1, 1, nx, ny)], ny * nz, MPI_DOUBLE, neighbors[4], 5, cart_comm, &requests[req_count++]);}
        if (neighbors[5] != MPI_PROC_NULL) { MPI_Irecv(&u[I3DMPI(nx + 1, 1, 1, nx, ny)], ny * nz, MPI_DOUBLE, neighbors[5], 4, cart_comm, &requests[req_count++]);}

        // SENDS
        if (neighbors[4] != MPI_PROC_NULL) { MPI_Isend(&u[I3DMPI(1, 1, 1, nx, ny)], ny * nz, MPI_DOUBLE, neighbors[4], 4, cart_comm, &requests[req_count++]);}
        if (neighbors[5] != MPI_PROC_NULL) { MPI_Isend(&u[I3DMPI(nx, 1, 1, nx, ny)], ny * nz, MPI_DOUBLE, neighbors[5], 5, cart_comm, &requests[req_count++]);}

        MPI_Waitall(req_count, requests, statuses);
    }

    // Update grid
    for (int i = 1; i <= nx; i++) 
    {
        for (int j = 1; j <= ny; j++) 
        {
            for (int k = 1; k <= nz; k++) 
            {
                u_new[I3DMPI(i, j, k, nx, ny)] = u[I3DMPI(i, j, k, nx, ny)]
                    + dt / dx2 * (u[I3DMPI(i + 1, j, k, nx, ny)] - 2 * u[I3DMPI(i, j, k, nx, ny)] + u[I3DMPI(i - 1, j, k, nx, ny)])
                    + dt / dy2 * (u[I3DMPI(i, j + 1, k, nx, ny)] - 2 * u[I3DMPI(i, j, k, nx, ny)] + u[I3DMPI(i, j - 1, k, nx, ny)])
                    + dt / dz2 * (u[I3DMPI(i, j, k + 1, nx, ny)] - 2 * u[I3DMPI(i, j, k, nx, ny)] + u[I3DMPI(i, j, k - 1, nx, ny)]);
            }
        }
    }
}


double compute_local_error_mpi(double *u, int nx, int ny, int nz, double T, double dx, double dy, double dz, int xo, int yo, int zo) 
{
    double local_error = 0.0;
    for (int i = 1; i <= nx; i++) 
    {
        for (int j = 1; j <= ny; j++) 
        {
            for (int k = 1; k <= nz; k++) 
            {
                double x = (i - 1 + xo) * dx;
                double y = (j - 1 + yo) * dy;
                double z = (k - 1 + zo) * dz;
                double exact = exp(-3 * T) * sin(x) * sin(y) * sin(z);
                double error = u[I3DMPI(i, j, k, nx, ny)] - exact;
                local_error += error * error * dx * dy * dz;
            }
        }
    }

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return local_error;
}


void initialize(double *u, int N, double dx, double dy, double dz) 
{
    for (int i = 0; i <= N; i++) 
    {
        for (int j = 0; j <= N; j++) 
        {
            for (int k = 0; k <= N; k++) 
            {
                double x = i * dx;
                double y = j * dy;
                double z = k * dz;
                u[I3D(i, j, k, N)] = sin(x) * sin(y) * sin(z);
            }
        }
    }
}


void time_step(double *u, double *u_new, int N, double dt, double dx2, double dy2, double dz2) 
{
    for (int i = 1; i < N; i++) 
    {
        for (int j = 1; j < N; j++) 
        {
            for (int k = 1; k < N; k++) 
            {
                u_new[I3D(i, j, k, N)] = u[I3D(i, j, k, N)] 
                    + dt / dx2 * (u[I3D(i+1, j, k, N)] - 2 * u[I3D(i, j, k, N)] + u[I3D(i-1, j, k, N)])
                    + dt / dy2 * (u[I3D(i, j+1, k, N)] - 2 * u[I3D(i, j, k, N)] + u[I3D(i, j-1, k, N)])
                    + dt / dz2 * (u[I3D(i, j, k+1, N)] - 2 * u[I3D(i, j, k, N)] + u[I3D(i, j, k-1, N)]);
            }
        }
    }
}


double compute_error(double *u, int N, double T, double dx, double dy, double dz) 
{
    double diff = 0.0;
    for (int i = 0; i <= N; i++) 
    {
        for (int j = 0; j <= N; j++) 
        {
            for (int k = 0; k <= N; k++) 
            {
                double x = i * dx;
                double y = j * dy;
                double z = k * dz;
                double exact = exp(-3 * T) * sin(x) * sin(y) * sin(z);
                double error = u[I3D(i, j, k, N)] - exact;
                diff += error * error * dx * dy * dz;
            }
        }
    }
    return sqrt(diff);
}
