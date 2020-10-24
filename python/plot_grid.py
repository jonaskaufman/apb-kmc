import sys
import numpy as np
from grid import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.optimize

def parse_profile_file(file_name):
    times = []
    profiles = [] 
    with open(file_name, 'r') as f:
        group = [] 
        while True:
            split_line = f.readline().split()
            if len(split_line) > 1:
                group.append(list(map(float, split_line)))
            else:
                if group:
                    profiles.append(group)
                if len(split_line) == 1:
                    times.append(float(split_line[0]))
                    group = []
                else:
                    break
    return times, profiles

# for grid values
def parse_grid_file(file_name):
    times = []
    grids = []
    with open(file_name, 'r') as f:
        next_grid = []
        while True:
            split_line = f.readline().split()
            if len(split_line) > 1:
                values = list(map(float, split_line))
                next_grid.append(values)
            else:
                if next_grid:
                    grids.append(
                        PixelGrid(len(next_grid[0]), len(next_grid)))
                    grids[-1].set_grid(next_grid)
                if len(split_line) == 1:
                    times.append(float(split_line[0]))
                    next_grid = []
                else:   # end of file
                    break
    return times, grids

phase_grid_file = sys.argv[1]
composition_grid_file = sys.argv[2]
composition_profile_file = sys.argv[3]
times, phase_grids = parse_grid_file(phase_grid_file)
times, composition_grids = parse_grid_file(composition_grid_file)
times, composition_profiles = parse_profile_file(composition_profile_file)
height = len(composition_profiles[0][0])
phase_height = 4*len(phase_grids[0].grid)
max_height = max(height, phase_height)
sigma = 40
smoothing_matrix = periodic_gaussian_matrix(height, sigma)

print(height/4)
print(np.mean(composition_profiles[0][0]))

for i in range(len(times)):
    fig, ax = plt.subplots(1 , 3, figsize=(4,7))
    ax[0].imshow(phase_grids[i].grid, aspect=2.0, origin='lower', vmin=-0.75, vmax=1.5, cmap='Greys')
    ax[1].imshow(composition_grids[i].grid, aspect=0.5, origin='lower', vmin=0, vmax=1, cmap='plasma')
    profile = composition_profiles[i][0]
    smooth_profile = smoothing_matrix.dot(profile)
    ax[0].axis('off')
    ax[1].axis('off')
    ax[0].set_ylim(0, max_height/4)
    ax[2].tick_params(direction='in', zorder=-1)
    ax[2].plot(profile, range(height), '0.9',zorder=0)
    ax[2].scatter(profile, range(height), s=1, c=profile, cmap='plasma', vmin=0, vmax=1, zorder=1, clip_on=False)
    ax[2].plot(smooth_profile, range(height), 'k')
    ax[1].set_ylim(0, max_height)
    ax[2].set_ylim(0, max_height)
#    ax[2].set_xlim(-0.05, 0.55)
#    ax[2].set_xlim(-0.05, 1.05)   
    ax[0].set_xlabel('phase')
    ax[1].set_xlabel('composition')
    ax[2].set_xlabel('composition (average)')
    ax[2].spines['left'].set_visible(False)
    ax[2].spines['right'].set_visible(False)
    ax[2].spines['top'].set_visible(False)
    ax[2].get_yaxis().set_visible(False)
    plt.tight_layout()
    plt.show()
#    plt.savefig('plot.png', dpi=600)

times, composition_profiles = parse_profile_file(composition_profile_file)

