#ifndef WRAPPER_H
#define WRAPPER_H

#include "simulation.hpp"
#include <vector>

/// Calculate initial boundary spacings that will produce the desired composition profile
std::vector<int> spacings_for_sinusoidal_composition(const BOUNDARY_TYPE& boundary_type,
                                                     const double& composition_average,
                                                     const double& composition_amplitude,
                                                     const int& target_height);

/// Wrapper to set up and perform single or multiple simulations and process their output
class SimulationWrapper
{
public:
    SimulationWrapper() = delete;

    /// Construct a wrapper with given simulation conditions
    SimulationWrapper(const BOUNDARY_TYPE& boundary_type,
                      const int& width,
                      const std::vector<int>& initial_spacings,
                      const double& temperature);

    /// Perform a single simulation for a given number of passes
    //  The phase grid, composition grid, and composition profile are printed out
    //  (with a time stamp) every print_interval passes, to their respective file streams
    void perform_single(const int& total_passes,
                        const int& print_interval,
                        std::ofstream& phase_grid_file_stream,
                        std::ofstream& composition_grid_file_stream,
                        std::ofstream& composition_profile_file_stream);

    /// Perform a set of simulations under the same conditions, each for a given numberr of passes
    //  The composition profile for each simulation is printed out (with a time stamp)
    //  every print_interval passes, to the file stream
    void perform_set(const int& total_simulations,
                     const int& total_passes,
                     const int& print_interval,
                     std::ofstream& composition_profile_file_stream);

private:
    /// Boundary type for simulations
    const BOUNDARY_TYPE boundary_type;

    /// Width of simulation grid
    const int width;

    /// Initial boundary spacings for simulations (this determines the grid height)
    const std::vector<int> initial_spacings;

    /// Temperature for simulations
    const double temperature;

    /// Set up a single simulation
    Simulation setup() const;

    /// Determine whether grid should be staggered based on boundary_type
    bool stagger_grid() const;
};

#endif
