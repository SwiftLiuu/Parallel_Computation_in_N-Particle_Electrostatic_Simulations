import matplotlib.pyplot as plt
import numpy as np


grid_size = 10  
cell_size = 1   
highlight_cell = (4, 4)  

# Generate random particle positions within the grid
np.random.seed(1)
particles_x = np.random.uniform(0, grid_size, 230)
particles_y = np.random.uniform(0, grid_size, 230)

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(0, grid_size)
ax.set_ylim(0, grid_size)

for i in range(grid_size + 1):
    ax.axhline(i, color='gray', linestyle='--', linewidth=0.5)
    ax.axvline(i, color='gray', linestyle='--', linewidth=0.5)

for dx in [-1, 0, 1]:
    for dy in [-1, 0, 1]:
        x, y = highlight_cell[0] + dx, highlight_cell[1] + dy
        rect = plt.Rectangle((x, y), cell_size, cell_size, color='lightblue', alpha=0.3)
        ax.add_patch(rect)

ax.scatter(particles_x, particles_y, color='black', s=15, label='Particles')

red_particle_x = highlight_cell[0] + 0.5
red_particle_y = highlight_cell[1] + 0.5
ax.scatter(red_particle_x, red_particle_y, color='red', s=50, label='Current Particle Under Calculation')

ax.set_xlabel("X Position")
ax.set_ylabel("Y Position")
ax.set_title("Mode 1: Particle Interaction Within Local Grids")
ax.legend(loc="upper right")
plt.grid(False)
plt.show()
