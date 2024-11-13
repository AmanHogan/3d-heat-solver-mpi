import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("../out/mpi_results.csv", header=None, names=["Runtime", "Error", "nprocs"])

df["Runtime"] = pd.to_numeric(df["Runtime"], errors='coerce')
df["Error"] = pd.to_numeric(df["Error"], errors='coerce')
df["nprocs"] = pd.to_numeric(df["nprocs"], errors='coerce')

df.dropna(inplace=True)
runtimes = df["Runtime"]
nprocs = df["nprocs"]

print(runtimes)
print(nprocs)


ideal_runtimes = runtimes / nprocs
print(ideal_runtimes)
scaling_efficiency = (ideal_runtimes / runtimes) * 100
print(scaling_efficiency)

print("Scaling Efficiency:")
for idx, row in df.iterrows():
    print(f"Processors: {row['nprocs']}, Actual Runtime: {row['Runtime']}s, "
          f"Ideal Runtime: {ideal_runtimes[idx]:.3f}s, Efficiency: {scaling_efficiency[idx]:.2f}%")

plt.figure(figsize=(10, 6))
plt.plot(nprocs, runtimes, 'o-', label='Actual Runtime')
plt.plot(nprocs, ideal_runtimes, 'x--', label='Ideal (Perfect) Scaling')
plt.xlabel('Number of Processes')
plt.ylabel('Runtime (seconds)')
plt.title('Strong Scaling Performance')
plt.legend()
plt.grid(True)
plt.xticks(nprocs)
plt.savefig('results.png')
