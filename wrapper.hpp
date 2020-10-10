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

    ///
    SimulationWrapper(const BOUNDARY_TYPE& boundary_type,
                      const int& width,
                      const std::vector<int>& initial_spacings,
                      const double& temperature);

    /// Set up a single simulation
    Simulation setup();

    /// Perform a single simulation
    void perform_single(const int& total_steps, const int& print_interval, std::ofstream& output_file_stream);

    /// Perform a set of simulations under the same conditions
    void perform_set(const int& total_simulations, const int& total_steps, const int& print_interval, std::ofstream& output_file_stream);

private:
    /// Boundary type for simulations
    const BOUNDARY_TYPE boundary_type;

    /// Width of simulation grid
    const int width;

    /// Initial boundary spacings for simulations
    const std::vector<int> initial_spacings;

    /// Temperature for simulations
    const double temperature;

    /// Determine whether grid should be staggered based on boundary_type
    bool stagger_grid() const;
};

#endif
