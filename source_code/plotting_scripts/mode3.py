import matplotlib.pyplot as plt
import numpy as np

grid_size = 10     
cell_size = 1      
num_particles = 100  
num_processes = 4    

np.random.seed(2)
particles_x = np.random.uniform(0, grid_size, num_particles)
particles_y = np.random.uniform(0, grid_size, num_particles)

colors = ['blue', 'green', 'orange', 'purple']
process_assignment = np.random.choice(num_processes, num_particles)

fig, axs = plt.subplots(2, 2, figsize=(12, 12))
fig.suptitle("Mode 3: Particles Assigned to Processes with Global Particle Copy", fontsize=16)

for i in range(num_processes):
    ax = axs[i // 2, i % 2]
    ax.set_xlim(0, grid_size)
    ax.set_ylim(0, grid_size)
    ax.set_title(f"Process {i+1}")

    for j in range(grid_size + 1):
        ax.axhline(j, color='gray', linestyle='--', linewidth=0.5)
        ax.axvline(j, color='gray', linestyle='--', linewidth=0.5)

    ax.scatter(particles_x, particles_y, color='black', s=25)

    process_particles_x = particles_x[process_assignment == i]
    process_particles_y = particles_y[process_assignment == i]
    ax.scatter(process_particles_x, process_particles_y, color=colors[i], s=25, label=f'Process {i+1} Particles')

    ax.set_xlabel("X Position")
    ax.set_ylabel("Y Position")
    ax.legend(loc="upper right")

plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
plt.show()
