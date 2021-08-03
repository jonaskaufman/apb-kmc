#!/usr/bin/python3
# Perform Fourier analysis on the average composition profiles obtained from a set of KMC runs

from grid import *
from parse import *
import glob
import os
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


def get_fundamental_magnitudes(profiles):
    """ Extract the fundamental fourier magnitude for each profile, normalized by that of the first profile """
    fundamental_magnitudes = []
    for profile in profiles:
        fourier_transform = np.fft.fft(profile)
        fundamental = np.abs(fourier_transform[1] + fourier_transform[-1])
        fundamental_magnitudes.append(fundamental)
    fundamental_magnitudes = np.array(fundamental_magnitudes)
    fundamental_magnitudes = fundamental_magnitudes / fundamental_magnitudes[0]
    return fundamental_magnitudes


def exponential_decay(t, a, tau):
    return a*np.exp(-t/tau)


def main():
    # Parse simulation output
    composition_profile_file = 'composition_profile.out'
    simulation_directory_prefix = 'sim.*'
    all_profile_files = glob.glob(os.path.join(
        simulation_directory_prefix, composition_profile_file))
    all_times = []
    all_profiles = []
    for profile_file in all_profile_files:
        times, profiles = parse_profile_file(profile_file)
        assert len(times) == len(profiles)
        all_times.append(times)
        all_profiles.append(profiles)
    for times in all_times:
        assert len(times) == len(all_times[0])
    height = len(all_profiles[0][0])/grid_y_scaling
    print(f'height = {height}')
    all_times = np.array(all_times)

    # Fourier transform each set of profiles
    all_fundamental_magnitudes = []
    for profiles in all_profiles:
        all_fundamental_magnitudes.append(get_fundamental_magnitudes(profiles))
    all_fundamental_magnitudes = np.array(all_fundamental_magnitudes)

    # Average and fit to exponential decay
    avg_times = np.average(all_times, axis=0)[1:]
    std_times = np.std(all_times, axis=0)[1:]
    avg_fundamental_magnitudes = np.average(
        all_fundamental_magnitudes, axis=0)[1:]
    std_fundamental_magnitudes = np.std(all_fundamental_magnitudes, axis=0)[1:]
    print(f'max std dev of times: {np.max(std_times)}')
    print(f'max std dev of magnitudes: {np.max(std_fundamental_magnitudes)}')

    popt, pcov = curve_fit(exponential_decay, avg_times,
                           avg_fundamental_magnitudes, p0=[1.0, 1000], sigma=std_fundamental_magnitudes)
    perr = np.sqrt(np.diag(pcov))
    print(f'tau = {popt[1]} +/- {perr[1]}')
    #tau = popt[1]
    #D = (height**2)/(4*(np.pi**2)*tau)
    #print(f'D = {D}')
    fit_fundamental_magnitudes = [
        exponential_decay(t, *popt) for t in avg_times]

    alph = 0.25
    for i in range(len(all_times)):
        plt.plot(all_times[i], all_fundamental_magnitudes[i],
                 'tab:orange', alpha=alph, zorder=1)
    plt.errorbar(avg_times, avg_fundamental_magnitudes, xerr=std_times,
                 yerr=std_fundamental_magnitudes, color='tab:gray', zorder=2)
    plt.plot(avg_times, fit_fundamental_magnitudes, 'tab:blue', zorder=3)
    plt.xlabel('time')
    plt.ylabel('fundamental Fourier magnitude')
#    plt.yscale('log')
    plt.show()


if __name__ == "__main__":
    main()
