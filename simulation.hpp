#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>
#include <random>
#include <utility>
#include <vector>

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
class SimulationCellGrid
{
public:
    SimulationCellGrid() = delete;

    /// Construct grid of a given width with regions of alternating phase
    //  of heights given by spacings.
    //  Note that spacings must be of even length to ensure periodicity.
    SimulationCellGrid(const int& width, const std::vector<int>& spacings, bool stagger);

    /// Grid dimensions
    const int width;
    const int height;

    /// Get phase of cell located at x,y
    int get_cell_phase(const int& x, const int& y) const;

    /// Flip phase of cell located at x,y
    void flip_cell_phase(const int& x, const int& y);

    /// TODO change to get_phase_pixel_grid
    /// Get pixel representation of grid
    //  If stagger is false, pixel grid is a normal grid of square cells
    //  If stagger is true, cells are doubled along the x direction and staggered along y
    std::vector<std::vector<int>> get_pixel_grid() const;

    std::vector<std::vector<int>> get_spacings_pixel_grid() const;

    std::vector<int> get_spacings_horizontal_pixel_averages() const;

    // do this later
    /// Print pixel grid to stream
    void print_pixel_grid(std::ostream& stream) const;

private:
    /// Whether grid should be staggered in pixel representation
    bool stagger;

    /// Phase values of grid
    std::vector<std::vector<int>> phase_grid;

    /// Set phase of cell located at x,y
    void set_cell_phase(const int& x, const int& y, const int& phase);
};

/// Kinetic events are vectors of indices representing sites to flip
using Event = std::vector<std::pair<int, int>>;

/// Possible boundary types for simulation (zeta minus or zeta plus)
enum class BOUNDARY_TYPE
{
    MINUS,
    PLUS
};

/// Monte Carlo simulation
class Simulation
{
public:
    Simulation() = delete;
    Simulation(const BOUNDARY_TYPE& boundary_type, const SimulationCellGrid& initial_grid, const int& temperature);

    /// Attempt one event, update configuration and time appropriately
    void step();

    /// Print pixel grid to stream
    void print_grid(std::ostream& stream) const;

    /// Get simulation time
    double get_time() { return time; }

private:
    SimulationCellGrid grid;

    /// Type of boundaries being simulated
    BOUNDARY_TYPE boundary_type;

    /// Simulation time, in units of reciprocal vibrational prefactor
    double time;

    /// Simulation temperature, in kelvin
    int temperature;

    /// List of all possible events that could occur over the course of a simulation
    std::vector<Event> event_list;

    /// Random number generator
    RandomGenerator random_generator;

    /// Populate the event list based on grid dimensions and boundary type
    void populate_event_list();

    /// Calculate the rate of an event given the current configuration
    double calculate_rate(const Event& event) const;

    /// Calculate the repulsion energy change due to an event
    double calculate_total_repulsion_energy_change(const Event& event) const;
   
    /// 
    double calculate_repulsion_energy_change(const int& x, const int& y) const;
};


#endif
