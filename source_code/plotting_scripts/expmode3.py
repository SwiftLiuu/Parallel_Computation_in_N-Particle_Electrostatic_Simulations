import matplotlib.pyplot as plt
import numpy as np

process_groups = [
    "1 proc", "1 proc", "1 proc",
    "2 procs", "2 procs", "2 procs",
    "4 procs", "4 procs", "4 procs"
]
threads = [
    "8 threads", "16 threads", "22 threads",
    "4 threads", "8 threads", "11 threads",
    "2 threads", "4 threads", "5 threads"
]
execution_times = [
    55.290,  # 1 proc, 8 threads
    50.471,  # 1 proc, 16 threads 
    48.249,  # 1 proc, 22 threads (optimal for 1 proc)
    96.264,  # 2 procs, 4 threads
    91.829,  # 2 procs, 8 threads
    90.854,  # 2 procs, 11 threads (optimal for 2 procs)
    194.124,   # 4 procs, 2 threads
    177.776,  # 4 procs, 4 threads 
    155.432   # 4 procs, 5 threads (optimal for 4 procs)
]

colors = {
    "1 proc": "skyblue",
    "2 procs": "lightgreen",
    "4 procs": "salmon"
}

optimal_indices = [2, 5, 8]  

plt.figure(figsize=(12, 8))

x = np.arange(len(process_groups))  
bar_width = 0.6  

bars = plt.bar(x, execution_times, color=[colors[group] for group in process_groups], width=bar_width)

for idx in optimal_indices:
    bars[idx].set_hatch('//')

plt.xlabel("Configuration (Processes, Threads)")
plt.ylabel("Execution Time (seconds)")
#plt.title("Execution Time for Different MPI Configurations (Cutoff Radius = 17500)")

plt.xticks(x, [f"{proc}, {thread}" for proc, thread in zip(process_groups, threads)], rotation=45, ha="right")

plt.yticks(fontsize=12)
plt.xticks(x, [f"{proc}, {thread}" for proc, thread in zip(process_groups, threads)], rotation=45, ha="right", fontsize=12)

plt.tight_layout()
plt.show()
