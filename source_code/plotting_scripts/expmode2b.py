import matplotlib.pyplot as plt
import numpy as np

threads = np.arange(8, 23)
execution_time_17500 = [
    37.2807, 35.8547, 32.2724, 32.0240, 31.4691, 28.3081, 
    32.0610, 28.1769, 28.0981, 28.1568, 28.5388, 28.2455, 
    31.3836, 32.0015, 32.9094
]
execution_time_18000 = [
    36.9282, 33.8573, 32.2985, 30.9723, 29.8001, 30.1844, 
    29.8798, 29.524, 29.5036, 31.4519, 31.481, 30.0445, 
    29.9204, 29.6652, 30.7813
]

min_time_17500_index = np.argmin(execution_time_17500)
min_time_18000_index = np.argmin(execution_time_18000)

fig, ax = plt.subplots(figsize=(12, 6))

bar_width = 0.35
index = np.arange(len(threads))

# Plot for cutoff = 17500
for i, time in enumerate(execution_time_17500):
    if i == min_time_17500_index:
        ax.bar(index[i], time, bar_width, color='blue', alpha=0.8, linewidth=2, edgecolor='black', hatch='//', label='Cutoff Radius = 17500 (Fastest Time)' if i == 0 else "")
    else:
        ax.bar(index[i], time, bar_width, color='blue', alpha=0.8, linewidth=1, edgecolor='black', label='Cutoff Radius = 17500' if i == 0 else "")

# Plot for cutoff = 18000
for i, time in enumerate(execution_time_18000):
    if i == min_time_18000_index:
        ax.bar(index[i] + bar_width, time, bar_width, color='green', alpha=0.7, linewidth=2, edgecolor='black', hatch='//', label='Cutoff Radius = 18000 (Fastest Time)' if i == 0 else "")
    else:
        ax.bar(index[i] + bar_width, time, bar_width, color='green', alpha=0.7, linewidth=1, edgecolor='black', label='Cutoff Radius = 18000' if i == 0 else "")

ax.set_xlabel('Number of Threads', fontsize=14)
ax.set_ylabel('Execution Time (seconds)', fontsize=14)
ax.set_xticks(index + bar_width / 2)
ax.set_xticklabels(threads, fontsize=12)
ax.tick_params(axis='y', labelsize=12)  # Increase y-axis tick labels size
ax.legend(loc="upper right", bbox_to_anchor=(1.15, 1))

plt.tight_layout()
plt.show()
