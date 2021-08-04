#!/usr/local/bin/python3
# Perform Fourier analysis on the average composition profiles obtained from a set of KMC runs

from grid import *
from parse import *
import glob
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import interp1d
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


def get_interpolated_samples(y, x, x_sample):
    """ Interpolate y values linearly and return sampled values at given x """
    assert min(x_sample) >= min(x)
    assert max(x_sample) <= max(x)
    linear_interpolation = interp1d(x, y, kind='linear')
    return linear_interpolation(x_sample)


def exponential_decay(t, a, tau_inv):
    return a*np.exp(-t*tau_inv)


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

    # Interpolate and average
    min_sample_time = max(all_times[:, 0])
    max_sample_time = min(all_times[:, -1])
    n_samples = len(all_times[0])
    sample_times = np.linspace(min_sample_time, max_sample_time, n_samples)[1:]
    all_interp_fundamental_magnitudes = [get_interpolated_samples(
        all_fundamental_magnitudes[i], all_times[i], sample_times) for i in range(len(all_times))]
    avg_fundamental_magnitudes = np.average(
        all_interp_fundamental_magnitudes, axis=0)
    std_fundamental_magnitudes = np.std(
        all_interp_fundamental_magnitudes, axis=0)

    cutoff_time_index = len(sample_times) - 1
    for i, t in enumerate(sample_times):
        if (avg_fundamental_magnitudes[i] < 2*std_fundamental_magnitudes[i]):
            print(
                f'Average magnitude is below two standard deviations at time {t}')
            cutoff_time_index = i
            break

    # Fit to exponential decay
    inv_tau_guess = 1
    for i, t in enumerate(sample_times):
        if (avg_fundamental_magnitudes[i] < 1/np.e):
            inv_tau_guess = 1 / t
            break
    popt, pcov = curve_fit(exponential_decay, sample_times[:cutoff_time_index],
                           avg_fundamental_magnitudes[:cutoff_time_index], p0=[1.0, inv_tau_guess], sigma=std_fundamental_magnitudes[:cutoff_time_index], absolute_sigma=True)
    perr = np.sqrt(np.diag(pcov))
    print(f'prefactor = {popt[0]} +/- {perr[0]}')
    print(f'1/tau = {popt[1]} +/- {perr[1]}')
    print(f'tau = {1/popt[1]} +/- {perr[1]/popt[1]/popt[1]} ')

    # Plot results
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(4, 6))

    fit_fundamental_magnitudes = [
        exponential_decay(t, *popt) for t in sample_times[:cutoff_time_index]]
    alph = 0.25
    for i in range(len(all_times)):
        for a in ax:
            a.plot(all_times[i], all_fundamental_magnitudes[i],
                   'tab:orange', alpha=alph, zorder=1)
    for a in ax:
        a.axvline(sample_times[cutoff_time_index], color='k', zorder=2)
        a.errorbar(sample_times, avg_fundamental_magnitudes,
                   yerr=std_fundamental_magnitudes, color='tab:gray', zorder=3)
        a.plot(sample_times[:cutoff_time_index],
               fit_fundamental_magnitudes, 'tab:blue', zorder=4)

    plt.xlabel('time')
    for a in ax:
        a.set_ylabel('fundamental Fourier magnitude')
    ax[0].set_ylim([0, 1])
    ax[1].set_yscale('log')
    ax[1].set_ylim([1e-3, 1.5])
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
