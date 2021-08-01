#ifndef CALCULATOR_H
#define CALCULATOR_H

#include "definitions.hpp"
#include "grid.hpp"
#include <memory>

/// Boltzmann constant, in eV/K
#define KB 8.617333262145e-5

/// Defect energies in eV
#define PLUS_KINK_DEFECT 0.2204

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
#define PLUS_REPULSION 0.0504

/**
 * Calculator of event rates in kinetic Monte Carlo simulation
 */
class EventRateCalculator
{
public:
    EventRateCalculator() = delete;
    EventRateCalculator(BOUNDARY_TYPE boundary_type,
                        double temperature,
                        const std::shared_ptr<const std::vector<Event>>& event_list_ptr,
                        const std::shared_ptr<const CellGrid>& grid_ptr);

    /// Returns the rate of an event given its ID
    double calculate_rate(ID event_id) const;

private:
    /// Event list
    std::shared_ptr<const std::vector<Event>> event_list_ptr;

    /// Simulation grid
    std::shared_ptr<const CellGrid> grid_ptr;

    /// Type of boundary being simulated (true for minus, false for plus)
    const bool boundary_minus;

    /// Simulation temperature, in kelvin
    const double temperature;

    /// Returns true if event is allowed
    bool is_valid_event(const Event& event) const;
    bool is_valid_event_zeta_minus(const Event& event) const;
    bool is_valid_event_zeta_plus(const Event& event) const;

    /// Returns true if cell is at a boundary
    bool at_valid_boundary(const Coordinates& coordinates) const;

    /// Returns true if event passes additional considerations for zeta plus events
    bool passes_additional_checks_zeta_plus(const Event& event) const;

    /// Returns the kinetic barrier for an event
    double calculate_barrier(const Event& event) const;

    /// Returns the base barrier height and endpoint energy change for an event
    std::pair<double, double> calculate_base_barrier_and_energy_change_zeta_minus(const Event& event) const;
    std::pair<double, double> calculate_base_barrier_and_energy_change_zeta_plus(const Event& event) const;

    /// Returns the sublattice (A or B) of the atom in given cell based on boundary type and phase
    SUBLATTICE get_sublattice_of_cell(const Coordinates& coordinates) const;

    /// Returns the total boundary repulsion energy change for an event
    double calculate_repulsion_energy_change(const Event& event) const;
};

#endif
