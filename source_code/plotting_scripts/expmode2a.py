import matplotlib.pyplot as plt
threads = list(range(1, 23))
times = [
    244.412, 123.571, 82.9421, 62.3182, 51.2583, 46.806, 41.6998, 36.9282,
    33.8573, 32.2985, 30.9723, 29.8001, 30.1844, 29.8798, 29.524, 29.5036,
    31.4519, 31.481, 30.0445, 29.9204, 29.6652, 30.7813
]
serial_time = 247.116

speedup = [serial_time / t for t in times]

fig, ax1 = plt.subplots()

ax1.plot(threads, times, 'o-', color='blue', label='Execution Time (s)')
ax1.set_xlabel('Number of Threads')
ax1.set_ylabel('Time Taken (seconds)', color='blue')
ax1.tick_params(axis='y', labelcolor='blue')

ax2 = ax1.twinx()
ax2.plot(threads, speedup, 's-', color='green', label='Speedup')
ax2.set_ylabel('Speedup Compared to Serial Time', color='green')
ax2.tick_params(axis='y', labelcolor='green')

ax1.legend(loc='center right', bbox_to_anchor=(0.95, 0.4))
ax2.legend(loc='center right', bbox_to_anchor=(0.95, 0.3))

plt.grid()
plt.tight_layout()
plt.show()
