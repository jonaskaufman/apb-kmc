import sys
import numpy as np
from grid import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.optimize

def parse_horizontal_results(file_name):
    times = []
    horizontal_values = [] 
    with open(file_name, 'r') as f:
        group = [] 
        while True:
            split_line = f.readline().split()
            if len(split_line) > 1:
                group.append(list(map(float, split_line)))
            else:
                if group:
                    horizontal_values.append(group)
                if len(split_line) == 1:
                    times.append(float(split_line[0]))
                    group = []
                else:
                    break
    return times, horizontal_values

# for grid values
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

