# Perform Fourier analysis on the average composition profiles obtained from a set of KMC runs
# Requires a set of directories sim.* each containing a composition_profile.out file

from grid import *
from parse import *
import glob
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit


def get_real_fundamental_amplitudes(profiles):
    """Extract the real part of the fundamental amplitude for each profile, corresponding to the cosine component"""
    amplitudes = []
    for profile in profiles:
        real_fourier_transform = np.fft.rfft(profile)
        fundamental = np.real(real_fourier_transform[1]) / len(profile)
        amplitudes.append(fundamental)
    amplitudes = np.array(amplitudes)
    return amplitudes


def get_interpolated_samples(y, x, x_sample):
    """Interpolate y values linearly and return sampled values at given x"""
    assert min(x_sample) >= min(x)
    assert max(x_sample) <= max(x)
    linear_interpolation = interp1d(x, y, kind="linear")
    return linear_interpolation(x_sample)


def exponential_decay(t, a, tau_inv):
    """Decaying exponential evaluated at t, with prefactor a and inverse relaxation time tau_inv"""
    return a * np.exp(-t * tau_inv)


def main():
    # Parse simulation output
    composition_profile_file = "composition_profile.out"
    simulation_directory_prefix = "sim.*"
    all_profile_files = glob.glob(
        os.path.join(simulation_directory_prefix, composition_profile_file)
    )
    print(f"{len(all_profile_files)} simulations found")
    all_times = []
    all_profiles = []
    for profile_file in all_profile_files:
        times, profiles = parse_profile_file(profile_file)
        assert len(times) == len(profiles)
        all_times.append(times)
        all_profiles.append(profiles)
    for times in all_times:
        assert len(times) == len(all_times[0])
    height = len(all_profiles[0][0]) / grid_y_scaling
    print(f"height = {height}")
    all_times = np.array(all_times)

    # Fourier transform each set of profiles
    all_amplitudes = []
    for profiles in all_profiles:
        all_amplitudes.append(get_real_fundamental_amplitudes(profiles))
    all_amplitudes = np.array(all_amplitudes)

    # Interpolate and average
    min_sample_time = max(all_times[:, 0])
    max_sample_time = min(all_times[:, -1])
    n_samples = len(all_times[0])
    sample_times = np.linspace(min_sample_time, max_sample_time, n_samples)[1:]
    all_interp_amplitudes = [
        get_interpolated_samples(all_amplitudes[i], all_times[i], sample_times)
        for i in range(len(all_times))
    ]
    avg_amplitudes = np.average(all_interp_amplitudes, axis=0)
    std_amplitudes = np.std(all_interp_amplitudes, axis=0)
    initial_amplitude = avg_amplitudes[0]

    # Flip amplitudes if initial is negative
    if initial_amplitude < 0:
        all_amplitudes = [-a for a in all_amplitudes]
        avg_amplitudes = -avg_amplitudes
        initial_amplitude = -initial_amplitude
    print(f"initial amplitude: {initial_amplitude}")

    fit = True
    if fit:
        cut_index = len(sample_times)
        print(f"samples in fit: {cut_index}")

        # Fit to exponential decay
        inv_tau_guess = 0
        for i, t in enumerate(sample_times):
            if avg_amplitudes[i] < 1 / np.e:
                inv_tau_guess = 1 / t
                break
        guess = [initial_amplitude, inv_tau_guess]
        popt, pcov = curve_fit(
            exponential_decay,
            sample_times[:cut_index],
            avg_amplitudes[:cut_index],
            p0=guess,
            bounds=(0, np.inf),
            sigma=std_amplitudes[:cut_index],
            absolute_sigma=True,
        )
        perr = np.sqrt(np.diag(pcov))
        print(f"prefactor = {popt[0]} +/- {perr[0]}")
        print(f"1/tau = {popt[1]} +/- {perr[1]}")
        print(f"tau = {1/popt[1]} +/- {perr[1]/popt[1]/popt[1]} ")
        fit_amplitudes = exponential_decay(sample_times[:cut_index], *popt)

    # Plot results
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(4, 6))

    alph = 0.25
    for i in range(len(all_times)):
        for a in ax:
            a.plot(all_times[i], all_amplitudes[i], "tab:orange", alpha=alph, zorder=1)
    for a in ax:
        a.errorbar(
            sample_times,
            avg_amplitudes,
            yerr=std_amplitudes,
            color="tab:gray",
            zorder=3,
        )
        if fit:
            a.plot(
                sample_times[:cut_index],
                fit_amplitudes,
                "tab:blue",
                linestyle="-",
                zorder=4,
            )

    plt.xlabel("time")
    for a in ax:
        a.set_ylabel("amplitude")
    ax[1].set_yscale("log")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
