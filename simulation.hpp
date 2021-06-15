#ifndef SIMULATION_H
#define SIMULATION_H

#include "calculator.hpp"
#include "definitions.hpp"
#include "grid.hpp"
#include <fstream>
#include <iostream>
#include <lotto/rejection.hpp>
#include <lotto/rejection_free.hpp>
#include <random>
#include <vector>

/**
 * Kinetic Monte Carlo simulation
 */
class Simulation
{
public:
    Simulation() = delete;

    /// Constructs a simulation for given boundary type with initial grid and temperature (in kelvin)
    Simulation(BOUNDARY_TYPE boundary_type, const SimulationCellGrid& initial_grid, double temperature);

    /// Attempts one event, updates the time and configuration accordingly
    void step();

    /// Runs a number of steps equal to the number of events considered
    void pass();

    /// Returns the simulation time
    double get_time() const { return time; }

    // TODO: Check dimensions of grid outputs
    //
    /// Returns a pixel grid of phase values (0 or 1). Width is doubled for zeta minus simulation.
    PixelGrid get_phase_pixel_grid() const;

    /// Returns a pixel grid of composition values. Height is doubled. Width is doubled for zeta minus simulation.
    PixelGrid get_composition_pixel_grid() const;

    /// Returns average composition profile obtained by averaging each row of the composition pixel grid.
    std::vector<double> get_average_composition_profile() const;

private:
    /// Simulation cell grid to store phase values
    std::shared_ptr<SimulationCellGrid> grid_ptr;

    /// Rate calculator
    std::shared_ptr<EventRateCalculator> rate_calculator_ptr;

    /// List of all possible events that can occur over the course of a simulation
    std::shared_ptr<const std::vector<Event>> event_list_ptr;

    /// Event selector
    // TODO: Change to rejection-free eventually
    std::unique_ptr<lotto::RejectionEventSelector<ID, EventRateCalculator>> selector_ptr;

    /// Type of boundary being simulated (true for minus, false for plus)
    const bool boundary_minus;

    /// Simulation time, in units of reciprocal vibrational prefactor
    double time;

    /// Simulation temperature, in kelvin
    const double temperature;

    /// Returns the rate of an event given the current configuration
    double calculate_rate(const Event& event) const;

    /// Generate event list
    std::vector<Event> generate_event_list() const;

    // Generate event ID lis
    std::vector<ID> generate_event_id_list() const;

    // Generate impact table
    std::map<ID, std::vector<ID>> generate_impact_table() const;
};

#endif
