# 3D Heat Equation Solver

## Project Overview
This project implements a **3D Heat Equation Solver** using both a **single MPI process** and **multiple MPI processes**. The solver is optimized for performance and can run on a cluster of compute nodes.

The project includes:
- `3dhs.c`: Solves the heat equation using a single MPI process.
- `3dmpi.c`: Solves the heat equation using multiple MPI processes with domain decomposition.
- `extract.py` and `visualize.py`: Python utilities to extract results and generate visualizations.
- A shell script (`run_all.sh`) to compile, distribute, and execute the programs on the cluster.

## Project Structure
```
project-root/
├── bin/ # Compiled executables
├── src/ # Source files (.c and .h)
├── include/ # Header files (.h)
├── scripts/ # Python and shell scripts
├── data/# Hosts file and other data
├── out/ # Output files
└── README.md # Project documentation
```

## Requirements
To compile and run this project, you need:
- `mpicc` (MPI C compiler)
- `mpirun` (MPI runtime)
- Python 3.x with `matplotlib`, `numpy`, and `pandas` for visualization

## Cluster Environment Setup
Before running the project on a cluster, you need to set up passwordless SSH access between nodes and ensure that all nodes are using Bash as the default shell. You can follow the instructions given by the profesor and TA.

## How to Run

### Setup script and hosts
For my run I used these hosts in the hostfile inside of the `scripts` folder:
```
192.168.5.132 slots=28
192.168.5.133 slots=28
```
Look at the top of the `run_all.sh script` and make sure the hosts in this file match the ones in the hostfile. These should be set correctly by default without you having to change anything.

Also of note, the bell cluster cannot really handle grid sizes of 150 >= in a reasonable amount of time. So the default values are set to lower than that value.

### Run the programs
First, do `cd scripts` in your terminal.

The project can be compiled and run using the provided `run_all.sh` script.

Enter into the console: `chmod 755 run_all.sh` then
`./run_all.sh`.

This will run the single version of the 3d heat solver, and the mpi version, and it will run a python script to give a csv output. 

## Running without script

### Running parallel version without script
If you just want to run the parallel version, without the script, you can do:

```bash
cd scripts
```
compile 

```bash
mpicc -O3 -march=native -funroll-loops -o ../bin/3dmpi ../src/3dmpi.c ../src/3dsolver.c -lm
```
then run

```bash
mpirun -np <no proccess> ../bin/3dmpi <x> <y> <z>
```

### Running sequential version without script

```bash
cd scripts
```
compile 

```bash
mpicc -O3 -march=native -funroll-loops -o ../bin/3dhs ../src/3dhs.c ../src/3dsolver.c -lm
```
then run

```bash
mpirun -np 1 ../bin/3dhs
```

## Output Files
The output from the sequential and parallel programs will be saved to:
- `out/out_seq.txt`
- `out/out_parallel.txt`

## Data Extraction and Visualization
After running the simulations, you can extract the results and generate visualizations using the provided Python scripts:

**Extract data**:
```bash
python3 scripts/extract.py
```

**Visualize the results**:
```bash
python3 scripts/visualize.py
```

## Authors
- Aman Hogan-Bailey
- The University of Texas at Arlington