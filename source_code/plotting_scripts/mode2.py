import matplotlib.pyplot as plt
import numpy as np

grid_size = 10  
cell_size = 1   
num_particles = 250  
num_threads = 4  

np.random.seed(2)
particles_x = np.random.uniform(0, grid_size, num_particles)
particles_y = np.random.uniform(0, grid_size, num_particles)

colors = ['blue', 'green', 'orange', 'purple']
thread_assignment = np.random.choice(num_threads, num_particles)

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(0, grid_size)
ax.set_ylim(0, grid_size)

for i in range(grid_size + 1):
    ax.axhline(i, color='gray', linestyle='--', linewidth=0.5)
    ax.axvline(i, color='gray', linestyle='--', linewidth=0.5)

for i in range(num_threads):
    thread_particles_x = particles_x[thread_assignment == i]
    thread_particles_y = particles_y[thread_assignment == i]
    ax.scatter(thread_particles_x, thread_particles_y, color=colors[i], s=25, label=f'Thread {i+1}')

ax.set_xlabel("X Position")
ax.set_ylabel("Y Position")
ax.set_title("Mode 2: Particle Calculation Assigned to Threads")
ax.legend(loc="upper right")
plt.grid(False)
plt.show()
