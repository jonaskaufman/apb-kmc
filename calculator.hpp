#ifndef CALCULATOR_H
#define CALCULATOR_H

#include "definitions.hpp"
#include "grid.hpp"

/// Boltzmann constant, in eV/K
#define KB 8.617333262145e-5

/// Hop barriers in eV
#define MINUS_KINK_FORM 0.0297
#define MINUS_KINK_MOVE 0.0327
#define PLUS_KINK_FORM_I 0.2766
#define PLUS_KINK_DESTROY_I 0.0563
#define PLUS_KINK_FORM_II 0.2342
#define PLUS_KINK_DESTROY_II 0.0139
#define PLUS_KINK_MOVE_I 0.0415
#define PLUS_KINK_MOVE_II 0.0532
#define PLUS_KINK_MOVE_III 0.0152
#define PLUS_KINK_MOVE_IV 0.1008

/// Nearest neighbor boundary repulsion energy in eV / boundary length unit
#define MINUS_REPULSION 0.0543
// TODO add repulsion!
#define PLUS_REPULSION 0.0

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

    ///
    bool at_valid_boundary(const std::pair<int, int>& coordinates, const SimulationCellGrid& grid) const;

    ///
    bool passes_additional_checks_zeta_plus(const Event& event, const SimulationCellGrid& grid) const;

    /// Returns the kinetic barrier for an event
    double calculate_barrier(const Event& event, const SimulationCellGrid& grid) const;

    ///
    double calculate_barrier_zeta_minus(const Event& event, const SimulationCellGrid& grid) const;

    ///
    double calculate_barrier_zeta_plus(const Event& event, const SimulationCellGrid& grid) const;

    /// Returns the sublattice (A or B) of the atom in given cell based on boundary type and phase
    SUBLATTICE get_sublattice_of_cell(const std::pair<int, int>& coordinates, const SimulationCellGrid& grid) const;

    /// Returns the total boundary repulsion energy change for an event
    double calculate_total_repulsion_energy_change(const Event& event, const SimulationCellGrid& grid) const;

    /// Returns the boundary repulsion energy change for a single cell flip
    double calculate_repulsion_energy_change(int x, int y, const SimulationCellGrid& grid) const;
};

#endif
