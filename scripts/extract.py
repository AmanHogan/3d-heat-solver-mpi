############################
# Author: Aman Hogan-Bailey
# Converts text output to csv file
# for proccessing
############################

import re
import csv

# Input and output file paths
input_file = '../out/out_parallel.txt'
output_csv = '../out/mpi_results.csv'

# Initialize lists to store data
runtimes = []
errors = []
nprocs_list = []
process_info = {}

# Regular expressions to extract data
process_count_pattern = re.compile(r"INFO_START: Using (\d+) proccesses")
process_info_pattern = re.compile(r"INFO: Process (\d+): \(x, y, z\) = \((\d+), (\d+), (\d+)\)")
error_pattern = re.compile(r"GLOBAL_ERROR: ([\deE\.\-]+)")
runtime_pattern = re.compile(r"RUNTIME: ([\d\.]+) seconds")

# Read the input file
with open(input_file, 'r') as file:
    data = file.read()

# Extract data
process_counts = process_count_pattern.findall(data)
errors = [float(e) for e in error_pattern.findall(data)]
runtimes = [float(r) for r in runtime_pattern.findall(data)]
process_matches = process_info_pattern.findall(data)

# Sort and deduplicate process counts
process_counts = sorted([int(x) for x in list(set(process_counts))])

# Extract process information
for match in process_matches:
    process_id, x, y, z = match
    process_info[int(process_id)] = (int(x), int(y), int(z))

# Ensure we have consistent data lengths
num_runs = min(len(process_counts), len(errors), len(runtimes))
process_counts = process_counts[:num_runs]
errors = errors[:num_runs]
runtimes = runtimes[:num_runs]

# Create csv
with open(output_csv, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Runtime (s)', 'Global Error', 'nprocs'])
    for i in range(num_runs):
        writer.writerow([runtimes[i], errors[i], process_counts[i]])
print(f"Data extracted and saved to {output_csv}")

