#ifndef SIMULATION_H
#define SIMULATION_H

#include "grid.hpp"
#include <fstream>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

/**
 * Random number generator
 *
 * Allows sampling of random integers and reals.
 */
class RandomGenerator
{
public:
    /// Default constructor, automatically seeds generator
    RandomGenerator();

    /// Returns a random integer between 0 and maximum_value (inclusive)
    int sample_integer_range(int maximum_value);

    /// Returns a random real from the unit interval
    double sample_unit_interval();

private:
    /// 64-bit Mersenne Twister generator
    std::mt19937_64 generator;

    /// Uniform real distribution on unit interval
    std::uniform_real_distribution<double> uniform_unit_interval_distribution;
};

// Kinetic event represented by indices of grid cells to flip
using Event = std::vector<std::pair<int, int>>;

/// Antiphase boundary types (zeta minus or zeta plus)
enum class BOUNDARY_TYPE
{
    MINUS,
    PLUS
};

/// Sublattices on honeycomb network
enum class SUBLATTICE
{
    A,
    B
};

/**
 * Kinetic Monte Carlo simulations
 */
class Simulation
{
public:
    Simulation() = delete;

    /// Constructs a simulation for given boundary type with initial grid and temperature (in kelvin)
    Simulation(BOUNDARY_TYPE boundary_type, const SimulationCellGrid& initial_grid, int temperature);

    /// Attempts one event, updates the time and configuration accordingly
    void step();

    /// Runs a number of steps equal to the number of events considered
    void pass();

    /// Returns a pixel grid of phase values (0 or 1). Width is doubled for zeta minus simulation.
    PixelGrid get_phase_pixel_grid() const;

    /// Returns a pixel grid of composition values. Height is doubled. Width is doubled for zeta minus simulation.
    PixelGrid get_composition_pixel_grid() const;

    /// Returns average composition profile obtained by averaging each row of the composition pixel grid.
    std::vector<double> get_average_composition_profile() const;

    /// Returns the simulation time
    double get_time() const { return time; }

private:
    /// Simulation cell grid to store phase values
    SimulationCellGrid grid;

    /// Type of boundary being simulated
    BOUNDARY_TYPE boundary_type;

    /// Simulation time, in units of reciprocal vibrational prefactor
    double time;

    /// Simulation temperature, in kelvin
    int temperature;

    /// List of all possible events that can occur over the course of a simulation
    std::vector<Event> event_list;

    /// Random number generator to be used throughout simulation
    RandomGenerator random_generator;

    /// Populates the event list based on the grid dimensions and boundary type
    void populate_event_list();

    /// Returns the sublattice (A or B) of the atom in given cell based on boundary type and phase
    SUBLATTICE get_sublattice_of_cell(int x, int y);

    /// Returns the rate of an event given the current configuration
    double calculate_rate(const Event& event) const;

    /// Returns the base kinetic barrier for an event
    double calculate_barrier(const Event& event) const;

    /// Returns the total boundary repulsion energy change for an event
    double calculate_total_repulsion_energy_change(const Event& event) const;

    /// Returns the boundary repulsion energy change for a single cell flip
    double calculate_repulsion_energy_change(int x, int y) const;
};

#endif
