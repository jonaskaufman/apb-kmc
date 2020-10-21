#ifndef CALCULATOR_H
#define CALCULATOR_H

#include "definitions.hpp"
#include "grid.hpp"

/// Boltzmann constant, in eV/K
#define KB 8.617333262145e-5

/// Barrier for hops in zeta minus ordering, in eV
//  TODO this is average, add in actual barriers
#define ZETA_MINUS_BARRIER 0.03

#define ZETA_MINUS_KINK_FORM 0.0297
#define ZETA_MINUS_KINK_MOVE_I 0.0324
#define ZETA_MINUS_KINK_MOVE_II 0.0331
#define ZETA_MINUS_REPULSION 0.0543

// TODO does this need to be a function?
double zeta_minus_boundary_energy(const int& spacing)
{
    switch (spacing)
    {
    case 1:
        return 0.0543;
    default:
        return 0;
    }
};


/**
 * Calculator of event rates in kinetic Monte Carlo simulation
 */
class EventRateCalculator
{
public:
    EventRateCalculator() = delete;
    EventRateCalculator(BOUNDARY_TYPE boundary_type, double temperature);

    /// Returns the rate of an event given the current configuration
    double calculate_rate(const Event& event, const SimulationCellGrid& grid) const;

private:
    BOUNDARY_TYPE boundary_type;

    double temperature;

    ///
    bool is_valid_event(const Event& event, const SimulationCellGrid& grid) const;

    ///
    bool is_valid_event_zeta_minus(const Event& event, const SimulationCellGrid& grid) const;

    ///
    bool is_valid_event_zeta_plus(const Event& event, const SimulationCellGrid& grid) const;

    /// Returns the kinetic barrier for an event
    double calculate_barrier(const Event& event, const SimulationCellGrid& grid) const;

    ///
    double calculate_barrier_zeta_minus(const Event& event, const SimulationCellGrid& grid) const;

    ///
    double calculate_barrier_zeta_plus(const Event& event, const SimulationCellGrid& grid) const;

    /// Returns the sublattice (A or B) of the atom in given cell based on boundary type and phase
    SUBLATTICE get_sublattice_of_cell(int x, int y, const SimulationCellGrid& grid);

    /// Returns the total boundary repulsion energy change for an event
    double calculate_total_repulsion_energy_change(const Event& event, const SimulationCellGrid& grid) const;

    /// Returns the boundary repulsion energy change for a single cell flip
    double calculate_repulsion_energy_change(int x, int y, const SimulationCellGrid& grid) const;
};

#endif
