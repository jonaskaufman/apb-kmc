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
sigma = 20
smoothing_matrix = periodic_gaussian_matrix(height, sigma)

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
    ax[2].set_xlim(-0.05, 0.55)   
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


"""
file_name = sys.argv[1]
times, horizontal_values = parse_horizontal_results(file_name)
n_simulations = len(horizontal_values[0])
height = len(horizontal_values[0][0])
true_height = height/2

sigma = 20
boundary_type = '-'
smoothing_matrix = periodic_gaussian_matrix(height, sigma)

colors = cm.plasma(np.linspace(0, 1, len(times)))
# get average profile for each frame
average_composition_profiles = []
for t in range(len(times)):
    profiles = []
    for k in range(n_simulations):
        # first gaussian smoothing
        smooth_profile = smoothing_matrix.dot(horizontal_values[t][k])
        if t == 0 and k == 0:
            plt.plot(horizontal_values[t][k])
            plt.plot(smooth_profile)
            plt.xlabel('position')
            plt.ylabel('composition')
            plt.show()
        profiles.append(smooth_profile)
    profiles = np.array(profiles)
    average_composition_profiles.append(np.mean(profiles, axis=0))
    plt.plot(average_composition_profiles[t], color=colors[t])
plt.ylim(0.4, 0.5)
plt.xlabel('position')
plt.ylabel('composition')
plt.show()

x_avg = np.mean(average_composition_profiles[0])
x_min = np.min(average_composition_profiles[0])
x_max = np.max(average_composition_profiles[0])
print(f'lambda = {true_height}')
print(f'x_avg = {x_avg}')
print(f'x_min = {x_min}')
print(f'x_max = {x_max}')
print(f'x_amp = {(x_max - x_min)/2}')
# TODO normalize FFT properly?
abs_ffts = []
for profile in average_composition_profiles:
    average_value = np.mean(profile)
    shifted_profile = profile - average_value
    fourier = np.fft.fft(shifted_profile)
    abs_ffts.append(np.abs(fourier))
abs_ffts = np.array(abs_ffts)
freq = np.fft.fftfreq(abs_ffts.shape[-1])

for i, time in enumerate(times):
    plt.plot(freq, abs_ffts[i], label=f'time = {times[i]}')
#plt.legend()
plt.show()

def exponential(x, a, k):
    return a*np.exp(x*k)

max_fft = np.amax(abs_ffts, axis=1)
supermax = np.amax(max_fft)
max_fft = max_fft/supermax
popt_exponential, pcov_exponential = scipy.optimize.curve_fit(exponential, times[1:], max_fft[1:], p0=[max(max_fft), -1/max(times)])
tau = -1/popt_exponential[1]
D = true_height*true_height/(4*np.pi*np.pi*tau) 
print(f'tau = {tau}')
print(f'D = {D}')
fit_values = [exponential(t, *popt_exponential) for t in times[1:]]
plt.plot(times[1:], fit_values, 'k-')
plt.plot(times, max_fft, 'o')
plt.yscale('log')
plt.show()
"""
