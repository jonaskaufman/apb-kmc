#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>
#include <fstream>
#include <random>
#include <utility>
#include <vector>

/// Random number generator
class RandomGenerator
{
public:
    /// Default constructor
    RandomGenerator();

    /// Get random integer between 0 and maximum_value (inclusive)
    int sample_integer_range(const int& maximum_value);

    /// Get a random real in the unit interval
    double sample_unit_interval();

private:
    /// 64-bit Mersenne Twister generator
    std::mt19937_64 generator;

    /// Uniform real distribution on unit interval
    std::uniform_real_distribution<double> uniform_unit_interval_distribution;
};

/// Periodic grid of simulation cells with phase 0 or 1
/// Can be staggered rectangular grid (like bricks) or square grid
class SimulationCellGrid
{
public:
    SimulationCellGrid() = delete;

    /// Construct grid of a given width with regions of alternating phase
    //  of heights given by spacings.
    //  An even number of spacings must be given to ensure periodicity
    SimulationCellGrid(const int& width, const std::vector<int>& spacings, bool stagger);

    /// Grid dimensions
    const int width;
    const int height;

    /// Get phase of cell located at x,y
    bool get_cell_phase(const int& x, const int& y) const;

    /// Flip phase of cell located at x,y from 0 to 1, or vice versa
    void flip_cell_phase(const int& x, const int& y);

    // TODO: Use typdef or struct for pixel grids?
    /// Get pixel representation of phase grid
    //  If stagger is false, pixel grid is a normal grid of square cells
    //  If stagger is true, cells are doubled along the x direction and staggered along y
    std::vector<std::vector<int>> get_phase_pixel_grid() const;

    /// Get pixel representation of grid where the value of each pixel
    //  is equal to the combined vertical distance to the nearest two boundaries
    std::vector<std::vector<int>> get_spacings_pixel_grid() const;

private:
    /// Whether grid should be staggered in pixel representation
    bool stagger;

    /// Phase values of grid
    std::vector<std::vector<bool>> phase_grid;

    /// Set phase of cell located at x,y
    void set_cell_phase(const int& x, const int& y, const bool& phase);
};

// TODO: Come up with a better name
/// Average pixel values across each row to obtain a vertical profile
std::vector<double> get_horizontal_pixel_averages(const std::vector<std::vector<int>>& pixel_grid);

/// Print pixel grid to output stream
void print_pixel_grid(const std::vector<std::vector<int>>& pixel_grid, std::ostream& stream);

// Kinetic events are vectors of indices representing grid cells to flip
using Event = std::vector<std::pair<int, int>>;

/// Possible boundary types for simulation (zeta minus or zeta plus)
enum class BOUNDARY_TYPE
{
    MINUS,
    PLUS
};

/// Kinetic Monte Carlo simulation
class Simulation
{
public:
    Simulation() = delete;
    Simulation(const BOUNDARY_TYPE& boundary_type, const SimulationCellGrid& initial_grid, const int& temperature);

    /// Attempt one event, update time, update configuration (if accepted)
    void step();

    /// Print phase pixel grid to stream
    void print_phase_pixel_grid(std::ostream& stream) const;

    /// Print spacings pixel grid to stream
    void print_spacings_pixel_grid(std::ostream& stream) const;

    ///
    std::vector<double> get_horizontal_pixel_average_spacings() const;

    /// Get simulation time
    double get_time() const { return time; }

private:
    /// Cell grid storing phase values
    SimulationCellGrid grid;

    /// Type of boundaries being simulated
    BOUNDARY_TYPE boundary_type;

    /// Simulation time, in units of reciprocal vibrational prefactor
    double time;

    /// Simulation temperature, in kelvin
    int temperature;

    /// List of all possible events that could occur over the course of a simulation
    std::vector<Event> event_list;

    /// Random number generator to be used throughout simulation
    RandomGenerator random_generator;

    /// Populate the event list based on grid dimensions and boundary type
    void populate_event_list();

    /// Calculate the rate of an event given the current configuration
    double calculate_rate(const Event& event) const;

    /// Calculate base barrier for event
    double calculate_barrier(const Event& event) const;

    /// Calculate the repulsion energy change due to an event
    double calculate_total_repulsion_energy_change(const Event& event) const;

    /// Calculate the repulsion energy change due to a single site flip
    double calculate_repulsion_energy_change(const int& x, const int& y) const;
};

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
