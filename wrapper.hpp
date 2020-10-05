#ifndef WRAPPER_H
#define WRAPPER_H

#include <vector>
#include "simulation.hpp"

// Calculate initial spacings from composition profile
std::vector<int> spacings_for_sinusoidal_composition(const BOUNDARY_TYPE& boundary_type,
                                                     const double& composition_average,
                                                     const double& composition_amplitude,
                                                     const int& target_height);

/// Wrapper to set up and perform single or multiple simulations and process their output
class SimulationWrapper
{
public:
    SimulationWrapper() = delete;
    SimulationWrapper(const BOUNDARY_TYPE& boundary_type,
                      const std::vector<int>& initial_spacings,
                      const double& temperature);
    void perform_single(const int& total_steps, const int& print_interval, const std::ofstream& output_file_stream);
    void perform_set(const int& total_steps, const int& print_interval, const std::ofstream& output_file_stream);

private:
    const BOUNDARY_TYPE boundary_type;
    const std::vector<int> initial_spacings;
    const double temperature;
    //    Simulation simulation;
};

#endif
