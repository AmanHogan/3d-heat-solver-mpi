/**
 * Author: Aman Hogan-Bailey
 * Header file that contains function definitions for 
 * 3d heatsolver for both mpi and mpi on one proccess
 */

#ifndef THREED_SOLVER_H
#define THREED_SOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

// Constants
#define PI 3.14159265358979323846

// Indexing for 3d heat solver using mpi, accounting for ghost cells
#define I3DMPI(i, j, k, nx, ny) ((i) + (j) * (nx + 2) + (k) * (nx + 2) * (ny + 2))

// Indexing for 3d heat solver
#define I3D(i, j, k, N) ((i) + (j) * (N + 1) + (k) * (N + 1) * (N + 1))

/**
 * Initializes the 3D grid with the condition: `u(x,y,z) = sin(x) * sin(y) * sin(z)` 
 * for use in MPI communication using ghost cells.
 * @param u 3D grid array
 * @param nx grid points in x
 * @param ny grid points in y
 * @param nz grid points in z
 * @param dx spacing x
 * @param dy spacing y
 * @param dz spacing z
 * @param xo offset in x
 * @param yo offset in y
 * @param zo offset in z
 * @returns void
 * 
 */
void initialize_mpi(double *u, int nx, int ny, int nz, double dx, double dy, double dz, int xo, int yo, int zo);

/**
 * Initializes the 3D grid with the condition: `u(x,y,z) = sin(x) * sin(y) * sin(z)` 
 * for use in MPI communication with only just `one` proccess
 * @param u 3D grid array
 * @param N num of grid points
 * @param dx spacing x
 * @param dy spacing y
 * @param dz spacing z
 * @returns void
 * 
 */
void initialize(double *u, int N, double dx, double dy, double dz);

/**
 * Performs one time step for the MPI-based heat solver.
 * @param u current D grid
 * @param u_new updated 3D grid after the time step.
 * @param nx grid points in x
 * @param ny grid points in y
 * @param nz grid points in z
 * @param dt time step
 * @param dx2 spacing x squared
 * @param dy2 spacing y squared
 * @param dz2 spacing z squared
 * @param cart_comm cartisean communicatior
 * @param neighbors vector containing ranks of differnt procccess in 3d directions
 * @returns void
 */
void time_step_mpi(double *u, double *u_new, int nx, int ny, int nz, double dt, double dx2, double dy2, double dz2, MPI_Comm cart_comm, int neighbors[6]);

/**
 * Performs one time step for the 3d heat solver.
 * @param u current D grid
 * @param u_new updated 3D grid after the time step.
 * @param N no of points
 * @param dt time step
 * @param dx2 spacing x squared
 * @param dy2 spacing y squared
 * @param dz2 spacing z squared
 * @returns void
 */
void time_step(double *u, double *u_new, int N, double dt, double dx2, double dy2, double dz2);

/**
 * Computes local error of a given proccess and returns that value.
 * @param u 3D grid array
 * @param nx grid points in x
 * @param ny grid points in y
 * @param nz grid points in z
 * @param T end time
 * @param dx spacing x
 * @param dy spacing y
 * @param dz spacing z
 * @param xo offset in x
 * @param yo offset in y
 * @param zo offset in z
 * @returns local_error
 */
double compute_local_error_mpi(double *u, int nx, int ny, int nz, double T, double dx, double dy, double dz, int xo, int yo, int zo);

/**
 * Computes local error and returns that value
 * @param u 3D grid array
 * @param N no of points
 * @param T end time
 * @param dx spacing x
 * @param dy spacing y
 * @param dz spacing z
 * @returns local_error
 */
double compute_error(double *u, int N, double T, double dx, double dy, double dz);

#endif // THREED_SOLVER_H
