# Plot phase grid, composition grid, and composition profile for each snapshot of a single KMC run
# Requires phase_grid.out, composition_grid.out, and composition_profile.out files

from grid import *
from parse import *
import matplotlib.pyplot as plt


def plot_and_report(
    time, phase_grid, composition_grid, composition_profile, smooth_composition_profile
):
    """Plot the phase/composition grids and original/smoothed composition profiles, report time and composition information"""

    # TODO make these calculations into common functions
    x_min = min(smooth_composition_profile)
    x_max = max(smooth_composition_profile)
    x_ampl = (x_max - x_min) / 2
    x_avg = np.mean(smooth_composition_profile)
    print(f"time = {time}\tx_avg = {x_avg}\t x_ampl = {x_ampl}")
    # Determine limits, aspect ratio
    phase_height = len(phase_grid)
    composition_height = len(composition_grid)
    max_height = max(phase_height, composition_height)

    if x_max > 0.5:
        composition_min = 0.5
        composition_max = 1.0
    else:
        composition_min = 0.0
        composition_max = 0.5
    composition_pad = 0.05
    aspect_ratio = grid_x_scaling / grid_y_scaling

    # Set up subplots
    fig, ax = plt.subplots(1, 3, figsize=(2.5, 3.5), sharey=True)

    pvmin = -0.25
    pvmax = 1.25
    pcmap = "viridis_r"
    # Plot grids
    ax[0].imshow(
        phase_grid,
        aspect=aspect_ratio,
        origin="lower",
        vmin=-0.75,
        vmax=2.00,
        cmap="Greys",
    )
    ax[1].imshow(
        composition_grid,
        aspect=aspect_ratio,
        origin="lower",
        vmin=pvmin,
        vmax=pvmax,
        cmap=pcmap,
    )
    ax[0].set_ylim(0, max_height)
    ax[0].set_xlabel("phase")
    ax[1].set_xlabel("composition")
    ax[0].axis("off")
    ax[1].axis("off")

    # Plot profiles
    heights = range(composition_height)
    ax[2].scatter(
        composition_profile,
        heights,
        s=1,
        c=composition_profile,
        cmap=pcmap,
        vmin=pvmin,
        vmax=pvmax,
        zorder=1,
        clip_on=False,
    )
    ax[2].plot(smooth_composition_profile, heights, "k")
    ax[2].set_xlim(composition_min - composition_pad, composition_max + composition_pad)
    ax[2].set_xlabel(r"$x$")
    ax[2].tick_params(direction="in", zorder=-1)
    ax[2].spines["left"].set_visible(False)
    ax[2].spines["right"].set_visible(False)
    ax[2].spines["top"].set_visible(False)
    ax[2].get_yaxis().set_visible(False)

    # Show
    plt.tight_layout()
    plt.show()


def main():
    # Parse simulation output
    phase_grid_file = "phase_grid.out"
    composition_grid_file = "composition_grid.out"
    composition_profile_file = "composition_profile.out"
    times, phase_grids = parse_grid_file(phase_grid_file)
    times, composition_grids = parse_grid_file(composition_grid_file)
    times, composition_profiles = parse_profile_file(composition_profile_file)
    assert len(times) == len(phase_grids)
    assert len(times) == len(composition_grids)
    assert len(times) == len(composition_profiles)
    height = len(composition_profiles[0]) / grid_y_scaling
    print(f"height = {height}")

    # Smooth composition profile
    smooth_composition_profiles = periodic_smooth_profiles(
        composition_profiles, y_sigma_scaled
    )

    # Plot and report for each time step
    for i in range(len(times)):
        plot_and_report(
            times[i],
            phase_grids[i],
            composition_grids[i],
            composition_profiles[i],
            smooth_composition_profiles[i],
        )


if __name__ == "__main__":
    main()
