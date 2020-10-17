#ifndef SIMULATION_H
#define SIMULATION_H

// TODO reduce includes
#include <fstream>
#include <iostream>
#include <random>
#include <utility>
#include <vector>
#include "grid.hpp"

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

// Kinetic events are vectors of indices representing grid cells to flip
using Event = std::vector<std::pair<int, int>>;

/// Possible boundary types for simulation (zeta minus or zeta plus)
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

/// Kinetic Monte Carlo simulation
class Simulation
{
public:
    Simulation() = delete;
    Simulation(const BOUNDARY_TYPE& boundary_type, const SimulationCellGrid& initial_grid, const int& temperature);

    /// Attempt one event, update time, update configuration (if accepted)
    void step();

    /// Attempt a number of steps equal to the total number of events considered 
    void pass();

    ///
    PixelGrid get_phase_pixel_grid() const;

    ///
    PixelGrid get_composition_pixel_grid() const;

    ///
    std::vector<double> average_horizontal_composition_pixels() const;

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

    /// Get the sublattice of atom in given cell, based on boundary type and phase
    SUBLATTICE get_sublattice_of_cell(const int& x, const int& y);
};

#endif
