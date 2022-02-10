# apb-kmc
2-D kinetic Monte Carlo (KMC) simulation of antiphase boundary (APB) motion.

## Background
This program simulates the macroscopic behavior resulting from a particular [cation diffusion mechanism](https://doi.org/10.1103/PhysRevMaterials.5.055401) identified in the layered material P3-Na<sub>*x*</sub>CoO<sub>2</sub>. In this mechanism, Na diffusion through the 2-D intercalation layers is facilitated by the motion of APBs. Two specific types of APBs, named zeta minus/plus, are supported.

For further details of the simulations and results, see [our publication](https://pubs.acs.org/doi/abs/10.1021/acs.chemmater.1c04152).

## Compilation
Be sure to clone the repository using `--recurse-submodules`.

This program relies on the [kmc-lotto](https://github.com/jonaskaufman/kmc-lotto) library, which is included as a submodule. Its repository contains [tests](https://github.com/jonaskaufman/kmc-lotto#installation) that you may wish to run.

Compile `apb-kmc` using the provided `Makefile` by running `make`. If you want to use a compiler other than `g++-9` you must edit the `Makefile`. The compiler must have C++17 support.

Internal simulation parameters (e.g. APB kink migration barriers) are hardcoded for P3-Na<sub>*x*</sub>CoO<sub>2</sub> but may be edited in `src/calculator.hpp`. Note that the program must be recompiled for any changes to take effect.

## Usage
To run a simulation, simply run the `apb-kmc` executable with no arguments. Settings are passed via an input file named `input.json`, which must be located in the working directory.

### Input
The `input.json` file must contain the following entries:
* `"boundary_type"`: `"-"` or `"+"` for zeta minus or zeta plus APBs, respectively.
* `"target_grid_dimensions"`: `[width, height]`, simulation grid dimensions in grid cell units. The grid width will equal the provided width, while the grid height will be close to the provided height.
* `"initialization"`
    * `"mode"`: Currently only `"sinusoidal"` is supported. This will initialize a roughly sinusoidal composition profile along the grid height, with wavelength equal to the grid height. 
    * `"composition_average"`: Target average *x* value of sinusoidal composition profile. Must be < 0.5 for zeta minus and > 0.5 for zeta plus.
    * `"composition_amplitude"`: Target amplitude of sinusoidal composition profile.
* `"temperature"`: Simulation temperature, in kelvin.
* `"total_passes"`: Total number of passes to run, where one pass involves running a number of KMC steps equal to the total number of possible events considered (this scales with the grid area).
* `"print_interval"`: Interval between simulation snapshots, in passes. 
* `"full_output"`: `true` or `false`, whether full grid information is written out at each snapshot, or just 1-D composition profiles (averaged along grid width).

### Output
Each output file contains the elapsed simulation time followed by some grid information, at each snapshot.
* `phase_grid.out`: 2-D grid of phase values (0 or 1), with the height added/removed by APBs excluded.
* `composition_grid.out`: 2-D grid of local composition values *x*, with the height added/removed by APBs included.
* `composition_profile.out`: 1-D composition profile (averaged along grid width), with the height added/removed by APBs included.

In the output files, each grid cell is represented by two horizontal pixels and four vertical pixels. Simulation times are given in units of inverse vibrational prefactor of the kinetic hops.

## Examples
Two examples with input and output files are provided.

### Single zeta plus simulation (`examples/zeta_plus_single`)
This is a single simulation of zeta plus APBs, with full output written. Running the provided `python/plot_run.py` script in this directory will plot the phase grid, composition grid, and composition profile at each snapshot.

### Set of zeta minus simulations (`examples/zeta_minus_set`)
This is a set of 10 independent simulations of zeta minus APBs, with only composition profiles written. Running the provided `python/analyze_run_set.py` script in this directory will plot the amplitude of each composition profile over time, along with the average amplitude, which is fit to a decaying exponential.
