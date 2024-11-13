#!/bin/bash

# Define source files for the sequential and parallel programs
SEQUENTIAL_SRC="../src/3dhs.c"
MPI_SRC="../src/3dmpi.c"
SOLVER_SRC="../src/3dsolver.c"
SOLVER_HEADER="../src/3dsolver.h"

SEQUENTIAL_EXE="../bin/3dhs"
MPI_EXE="../bin/3dmpi"
OUTPUT_SEQ="../out/out_seq.txt"
OUTPUT_MPI="../out/out_parallel.txt"
MPI_EXE_HOST="3dmpi"

# Define host IPs manually
HOST1="192.168.5.132"
HOST2="192.168.5.133"
HOSTFILE="hosts"

# Define MPI configurations
configs=("1 1 1" "2 2 2" "2 5 5")  # Different grid configurations
process_counts=(1 8 50)            # Corresponding process counts

# Clean up old files
rm -f $SEQUENTIAL_EXE $MPI_EXE $OUTPUT_SEQ $OUTPUT_MPI

echo "-------------------------------------------------------"
echo "Compiling the sequential program ($SEQUENTIAL_SRC)..."
mpicc -O3 -march=native -funroll-loops -o $SEQUENTIAL_EXE $SEQUENTIAL_SRC $SOLVER_SRC -lm
if [ $? -ne 0 ]; then
    echo "Error: Compilation of $SEQUENTIAL_SRC failed."
    exit 1
fi

echo "Compiling the parallel MPI program ($MPI_SRC)..."
mpicc -O3 -march=native -funroll-loops -o $MPI_EXE $MPI_SRC $SOLVER_SRC -lm
if [ $? -ne 0 ]; then
    echo "Error: Compilation of $MPI_SRC failed."
    exit 1
fi

# Ensure executables are executable
chmod +x ./$MPI_EXE

echo "-------------------------------------------------------"
echo "Running the sequential program..."
mpirun -np 1 ./$SEQUENTIAL_EXE > $OUTPUT_SEQ
echo "Sequential program finished. Output saved to $OUTPUT_SEQ."
echo "-------------------------------------------------------"

# Function to copy the executable to specified nodes and set permissions
copy_to_nodes() {
    local program=$1
    local remote_path=$(basename "$program")

    echo "Copying $program to nodes $HOST1 and $HOST2..."

    # Copy to HOST1
    scp "$program" "$HOST1:~/$remote_path"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to copy $program to $HOST1"
    else
        ssh "$HOST1" "chmod +x ~/$remote_path"
    fi

    # Copy to HOST2
    scp "$program" "$HOST2:~/$remote_path"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to copy $program to $HOST2"
    else
        ssh "$HOST2" "chmod +x ~/$remote_path"
    fi
}


echo "-------------------------------------------------------"
echo "Running the parallel MPI program..."

# Run each MPI configuration
for i in "${!configs[@]}"; do
    config="${configs[i]}"
    process_count="${process_counts[i]}"
    
    echo "-------------------------------------------------------"
    echo "Running configuration: MPIX=${config}"
    echo "Using ${process_count} processes"
    
    # Copy executable to all nodes before each configuration
    copy_to_nodes "$MPI_EXE"
    
    if [ "$process_count" -eq 50 ]; then
        # Use hostfile for the 50-process configuration
        mpirun --tag-output --hostfile $HOSTFILE -n "$process_count" $MPI_EXE_HOST $config >> $OUTPUT_MPI
    else
        mpirun -n "$process_count" $MPI_EXE $config >> $OUTPUT_MPI
    fi
done

echo "Parallel MPI program finished. Output saved to $OUTPUT_MPI."
echo "-------------------------------------------------------"

echo "Running Python extraction script..."
python3 extract.py
echo "Extraction complete. Results saved to mpi_results.csv"
echo "-------------------------------------------------------"
