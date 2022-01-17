#ifndef SIMULATION_H
#define SIMULATION_H

#include "calculator.hpp"
#include "definitions.hpp"
#include "grid.hpp"
#include <fstream>
#include <iostream>
#include <lotto/rejection_free.hpp>
#include <random>
#include <set>
#include <vector>

/**
 * Kinetic Monte Carlo simulation
 */
class Simulation
{
public:
    Simulation() = delete;

    /// Constructs a simulation for given boundary type with initial grid and temperature (in kelvin)
    Simulation(BOUNDARY_TYPE boundary_type, const CellGrid& initial_grid, double temperature);

    /// Attempts one event, updates the time and configuration accordingly
    void step();

    /// Runs a number of steps equal to the number of events considered
    void pass();

    /// Returns the simulation time
    double get_time() const { return time; }

    /// Returns a pixel grid of phase values (0 or 1).
    //  Cell width and height are scaled by factors of 2 and 4, respectively.
    PixelGrid get_phase_pixel_grid() const;

    /// Returns a pixel grid of composition values.
    //  Cell width and height are scaled by factors of 2 and 4, respectively.
    PixelGrid get_composition_pixel_grid() const;

    /// Returns average composition profile obtained by averaging each row of the composition pixel grid
    std::vector<double> get_average_composition_profile() const;

private:
    /// Simulation cell grid to store phase values
    std::shared_ptr<CellGrid> grid_ptr;

    /// Rate calculator
    std::shared_ptr<EventRateCalculator> rate_calculator_ptr;

    /// List of all possible events that can occur over the course of a simulation
    std::shared_ptr<const std::vector<Event>> event_list_ptr;

    /// Event selector
    std::unique_ptr<lotto::RejectionFreeEventSelector<ID, EventRateCalculator>> selector_ptr;

    /// Type of boundary being simulated (true for minus, false for plus)
    const bool boundary_minus;

    /// Simulation time, in units of reciprocal vibrational prefactor
    double time;

    /// Simulation temperature, in kelvin
    const double temperature;

    /// Returns the rate of an event given the current configuration
    double calculate_rate(const Event& event) const;

    /// Generates event list
    std::vector<Event> generate_event_list() const;

    /// Generates event ID list
    std::vector<ID> generate_event_id_list() const;

    /// Generates map from cell coordinates to IDs of events involving that cell
    std::map<Coordinates, std::vector<ID>> generate_coordinates_ids_map() const;

    /// Generates the list of cell coordinates that impact the rate of a given event
    std::set<Coordinates> generate_impact_neighborhood(const Event& impacted_event) const;

    /// Generates impact table
    std::map<ID, std::vector<ID>> generate_impact_table() const;
};

#endif
