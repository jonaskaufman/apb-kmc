#!/usr/bin/python3
# Perform Fourier analysis on the average composition profiles obtained from a set of KMC runs

from grid import *
from parse import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.optimize import curve_fit

def plot_profiles(profiles):
    colors = cm.plasma(np.linspace(0, 1, len(profiles)))
    max_height = len(profiles[0])/grid_y_scaling
    heights = np.linspace(0, max_height, len(profiles[0]))
    for i, p in enumerate(profiles):
        plt.plot(heights, p, color=colors[i])
    plt.show()

def exponential(x, a, k):
    return a*np.exp(x*k)

def main():
    # Parse simulation output
    composition_profile_file = "composition_profile.out"
    times, composition_profile_sets = parse_profile_file(composition_profile_file)
    assert len(times) == len(composition_profile_sets)

    # Average and smooth at each time
    average_composition_profiles = [np.mean(p_set, axis=0) for p_set in composition_profile_sets]
    smooth_average_composition_profiles = periodic_smooth_profiles(average_composition_profiles, y_sigma_scaled)

    height = len(smooth_average_composition_profiles[0])/grid_y_scaling
    print(f'height = {height}')

    # Plot profiles over time
    plot_profiles(smooth_average_composition_profiles)

    fft_magnitudes = []
    for profile in smooth_average_composition_profiles:
        shifted_profile = profile - np.mean(profile)
        fourier_transform = np.fft.fft(shifted_profile)
        fft_magnitudes.append(np.abs(fourier_transform))
    fft_magnitudes = np.array(fft_magnitudes)
    fft_magnitudes = fft_magnitudes/np.amax(fft_magnitudes) # normalize by maximum
    #fft_frequencies = np.fft.fftfreq(fft_magnitudes.shape[-1])
    #fft_frequencies = fft_frequencies*grid_y_scaling 
    #plt.plot(fft_frequencies, fft_magnitudes[0]) 
    #plt.show()

    fundamental_magntitude = np.amax(fft_magnitudes, axis=1)
    popt, pcov = curve_fit(exponential, times, fundamental_magntitude, p0=[1.0, 0.0])
    tau = -1/popt[1]
    D = (height**2)/(4*(np.pi**2)*tau) 
    print(f'tau = {tau}')
    print(f'D = {D}')
    fit_values = [exponential(t, *popt) for t in times]
    plt.plot(times, fit_values, 'k-')
    plt.plot(times, fundamental_magntitude, 'o')
    plt.yscale('log')
    plt.show()


if __name__ == "__main__":
    main()
