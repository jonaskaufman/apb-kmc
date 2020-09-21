import sys
import numpy as np
from grid import *
import matplotlib.pyplot as plt


def parse_results(file_name):
    times = []
    phase_grids = []
    with open(file_name, 'r') as f:
        next_grid = []
        while True:
            split_line = f.readline().split()
            if len(split_line) > 1:
                values = list(map(int, split_line))
                next_grid.append(values)
            else:
                if next_grid:
                    phase_grids.append(
                        PixelGrid(len(next_grid[0]), len(next_grid)))
                    phase_grids[-1].set_grid(next_grid)
                if len(split_line) == 1:
                    times.append(float(split_line[0]))
                    next_grid = []
                else:   # end of file
                    break
    return times, phase_grids


sigma = 10
boundary_type = '-'

n_runs = 10
profiles = []
for n in range(n_runs):
    print(n)
    file_name = f'../output/run_{n}.out'
    times, phase_grids = parse_results(file_name)
    profiles.append([])
    for i, phase_grid in enumerate(phase_grids):
        spacing_grid = get_spacing_pixel_grid(phase_grid)
        smooth_spacing_grid = spacing_grid.get_smooth_grid(sigma)
        composition_grid = get_composition_pixel_grid(
        smooth_spacing_grid, boundary_type)
        composition_profile = composition_grid.get_horizontal_averages()
        profiles[-1].append(composition_profile)
profiles = np.array(profiles)
avg_profiles = np.mean(profiles, axis=0)

abs_ffts = []
for profile in avg_profiles:
    shifted_profile = profile - np.mean(profile) 
    fourier = np.fft.fft(shifted_profile)
    abs_ffts.append(np.abs(fourier))
abs_ffts = np.array(abs_ffts)
print(abs_ffts.shape)
freq = np.fft.fftfreq(abs_ffts.shape[-1])

for i, time in enumerate(times):
    plt.plot(freq, abs_ffts[i], label=f'time = {times[i]}')

plt.legend()
plt.show()

max_fft = np.amax(abs_ffts, axis=1)
plt.plot(times, max_fft, 'o')
plt.yscale('log')
plt.show()
