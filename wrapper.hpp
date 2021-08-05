#ifndef WRAPPER_H
#define WRAPPER_H

#include "simulation.hpp"
#include <vector>

/// Calculate initial boundary spacings that will produce a roughly sinusoidal composition profile with close to the
// desired average, amplitude, and height (wavelength)
std::vector<int> spacings_for_sinusoidal_composition(BOUNDARY_TYPE boundary_type,
                                                     double composition_average,
                                                     double composition_amplitude,
                                                     double target_physical_height);

/// Average boundary spacing for a given boundary type and composition
double average_spacing_from_composition(BOUNDARY_TYPE boundary_type, double composition);

/// Average composition for a given boundary type and boundary spacing
double average_composition_from_spacing(BOUNDARY_TYPE boundary_type, double spacing);

/// Height accounting for height of boundaries themselves, given a set of boundary spacings
double physical_height(BOUNDARY_TYPE boundary_type, const std::vector<int>& boundary_cell_spacings);

/// Print a report summarizing the initial spacings
void print_initialization_report(BOUNDARY_TYPE boundary_type,
                                 const std::vector<int>& initial_spacings,
                                 double target_physical_height,
                                 double target_composition,
                                 std::ostream& output_stream);

/// Wrapper to set up and perform simulations and process their output
class SimulationWrapper
{
public:
    SimulationWrapper() = delete;

    /// Construct a wrapper with given simulation conditions
    SimulationWrapper(BOUNDARY_TYPE boundary_type,
                      int width,
                      const std::vector<int>& initial_spacings,
                      double temperature);

    /// Perform a single simulation for a given number of passes
    //  The phase grid, composition grid, and composition profile are printed out
    //  (with a time stamp) every print_interval passes, to their respective file streams
    //  If full_output is false, only composition profile is written
    void perform_single(int total_passes,
                        int print_interval,
                        bool full_output,
                        std::ofstream& phase_grid_file_stream,
                        std::ofstream& composition_grid_file_stream,
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
