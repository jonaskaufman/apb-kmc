#include "calculator.hpp"
#include <cmath>
#include <stdexcept>

inline double boltzmann_factor(double barrier, double temperature) { return std::exp(-barrier / (KB * temperature)); }

EventRateCalculator::EventRateCalculator(BOUNDARY_TYPE boundary_type,
                                         double temperature,
                                         const std::shared_ptr<const std::vector<Event>>& event_list_ptr,
                                         const std::shared_ptr<const SimulationCellGrid>& grid_ptr)
    : boundary_minus{boundary_type == BOUNDARY_TYPE::MINUS}, temperature{temperature}, event_list_ptr{event_list_ptr}, grid_ptr{grid_ptr}
{
}

double EventRateCalculator::calculate_rate(ID event_id) const
{
    const Event& event = event_list_ptr->at(event_id);
    if (!is_valid_event(event))
    {
        return 0.0;
    }
    else
    {
        double barrier = calculate_barrier(event);
        return boltzmann_factor(barrier, temperature);
    }
}

bool EventRateCalculator::is_valid_event(const Event& event) const
{
    if (boundary_minus)
    {
        return is_valid_event_zeta_minus(event);
    }
    else
    {
        return is_valid_event_zeta_plus(event);
    }
}

bool EventRateCalculator::is_valid_event_zeta_minus(const Event& event) const
{
    // Only single cell events are allowed
    if (event.size() != 1)
    {
        return false;
    }
    // Cell must be at a flat boundary and not be the only cell separating two boundaries
    return at_valid_boundary(event[0]);
}

bool EventRateCalculator::is_valid_event_zeta_plus(const Event& event) const
{
    // Only single or double cell events are allowed
    if ((event.size() != 1) && (event.size() != 2))
    {
        return false;
    }
    // Cell(s) must be at a flat boundary and not be the only cell separating two boundaries
    for (const auto& coordinates : event)
    {
        if (!at_valid_boundary(coordinates))
        {
            return false;
        }
    }
    // Additional checks, mostly ensuring no single cell kinks are introduced
    if (!passes_additional_checks_zeta_plus(event))
    {
        return false;
    }
    return true;
}

bool EventRateCalculator::passes_additional_checks_zeta_plus(const Event& event) const
{
    if (event.size() == 1)
    {
        bool phase = grid_ptr->get_cell_phase(event[0]);
        bool phase_W = grid_ptr->get_neighbor_phase(event[0], DIRECTION::W);
        bool phase_E = grid_ptr->get_neighbor_phase(event[0], DIRECTION::E);
        // Neighbors must be opposite phase from each other
        if (phase_W == phase_E)
        {
            // No single cell kinks should ever be present
            if (phase != phase_W)
            {
                throw std::runtime_error("Single atom kink detected in a zeta plus boundary, which is not allowed");
            }
            return false;
        }
        // No single cell kinks should be created
        else
        {
            bool phase_X = (phase == phase_W) ? phase_W : phase_E;
            DIRECTION direction_X = (phase == phase_W) ? DIRECTION::W : DIRECTION::E;
            bool phase_X_X = grid_ptr->get_neighbor_phase(grid_ptr->get_neighbor_indices(event[0], direction_X), direction_X);
            if (phase_X_X != phase_X)
            {
                return false;
            }
        }
    }
    else if (event.size() == 2)
    {
        // Event cells must be in expected order
        if (grid_ptr->get_neighbor_indices(event[0], DIRECTION::E) != event[1])
        {
            throw std::runtime_error("Cells for double atom event are in wrong order or not horizontally adjacent.");
        }
        bool phase_0 = grid_ptr->get_cell_phase(event[0]);
        bool phase_1 = grid_ptr->get_cell_phase(event[1]);
        // Event cells must have same phase
        if (phase_0 != phase_1)
        {
            return false;
        }
        // Neighbors must have same phase
        bool phase_0_W = grid_ptr->get_neighbor_phase(event[0], DIRECTION::W);
        bool phase_1_E = grid_ptr->get_neighbor_phase(event[1], DIRECTION::E);
        if (phase_0_W != phase_1_E)
        {
            return false;
        }
        // No single cell kinks should be created
        if (phase_0 == phase_0_W)
        {
            bool phase_0_W_W = grid_ptr->get_neighbor_phase(grid_ptr->get_neighbor_indices(event[0], DIRECTION::W), DIRECTION::W);
            bool phase_1_E_E = grid_ptr->get_neighbor_phase(grid_ptr->get_neighbor_indices(event[1], DIRECTION::E), DIRECTION::E);
            if ((phase_0_W_W != phase_0_W) || (phase_1_E_E != phase_1_E))
            {
                return false;
            }
        }
    }
    else
    {
        return false;
    }
    return true;
}

double EventRateCalculator::calculate_barrier(const Event& event) const
{
    auto barrier_and_energy_change =
        boundary_minus ? calculate_barrier_and_energy_change_zeta_minus(event) : calculate_barrier_and_energy_change_zeta_plus(event);
    double base_barrier = barrier_and_energy_change.first;
    double base_energy_change = barrier_and_energy_change.second;
    double repulsion_energy_change = calculate_repulsion_energy_change(event);
    double total_energy_change = base_energy_change + repulsion_energy_change;
    double barrier = (total_energy_change < 0) ? base_barrier : base_barrier + repulsion_energy_change;
    if (barrier < 0.0)
    {
        throw std::runtime_error("Barrier is negative.");
    }
    return barrier;
}

std::pair<double, double> EventRateCalculator::calculate_barrier_and_energy_change_zeta_minus(const Event& event) const
{
    double barrier = 0.0;
    double energy_change = 0.0;
    bool phase_W = grid_ptr->get_neighbor_phase(event[0], DIRECTION::W);
    bool phase_E = grid_ptr->get_neighbor_phase(event[0], DIRECTION::E);
    if (phase_W == phase_E)
    {
        barrier = MINUS_KINK_FORM;
    }
    else
    {
        barrier = MINUS_KINK_MOVE;
    }
    return std::make_pair(barrier, energy_change);
}

std::pair<double, double> EventRateCalculator::calculate_barrier_and_energy_change_zeta_plus(const Event& event) const
{
    double barrier = 0.0;
    double energy_change = 0.0;
    if (event.size() == 2)
    {
        bool phase = grid_ptr->get_cell_phase(event[0]);
        bool phase_sides = grid_ptr->get_neighbor_phase(event[0], DIRECTION::W);
        SUBLATTICE sublattice_0 = get_sublattice_of_cell(event[0]);
        if (phase == phase_sides)
        {
            energy_change = PLUS_KINK_DEFECT;
            barrier = ((sublattice_0 == SUBLATTICE::B) ? PLUS_KINK_FORM_I : PLUS_KINK_FORM_II);
        }
        else
        {
            energy_change = -PLUS_KINK_DEFECT;
            barrier = ((sublattice_0 == SUBLATTICE::A) ? PLUS_KINK_DESTROY_I : PLUS_KINK_DESTROY_II);
        }
    }
    else
    {
        bool phase_S = grid_ptr->get_neighbor_phase(event[0], DIRECTION::S);
        bool phase_E = grid_ptr->get_neighbor_phase(event[0], DIRECTION::E);
        SUBLATTICE sublattice_E = get_sublattice_of_cell(grid_ptr->get_neighbor_indices(event[0], DIRECTION::E));
        if (phase_S == phase_E)
        {
            barrier = ((sublattice_E == SUBLATTICE::A) ? PLUS_KINK_MOVE_I : PLUS_KINK_MOVE_II);
        }
        else
        {
            barrier = ((sublattice_E == SUBLATTICE::A) ? PLUS_KINK_MOVE_III : PLUS_KINK_MOVE_IV);
        }
    }
    return std::make_pair(barrier, energy_change);
}

bool EventRateCalculator::at_valid_boundary(const Coordinates& coordinates) const

{
    if (boundary_minus)
    {
        bool phase_NW = grid_ptr->get_neighbor_phase(coordinates, DIRECTION::NW);
        bool phase_NE = grid_ptr->get_neighbor_phase(coordinates, DIRECTION::NE);
        bool phase_SW = grid_ptr->get_neighbor_phase(coordinates, DIRECTION::SW);
        bool phase_SE = grid_ptr->get_neighbor_phase(coordinates, DIRECTION::SE);
        return ((phase_NW == phase_NE) && (phase_SW == phase_SE) && (phase_NW != phase_SW));
    }
    else
    {
        bool phase_N = grid_ptr->get_neighbor_phase(coordinates, DIRECTION::N);
        bool phase_NW = grid_ptr->get_neighbor_phase(coordinates, DIRECTION::NW);
        bool phase_NE = grid_ptr->get_neighbor_phase(coordinates, DIRECTION::NE);
        bool phase_S = grid_ptr->get_neighbor_phase(coordinates, DIRECTION::S);
        bool phase_SW = grid_ptr->get_neighbor_phase(coordinates, DIRECTION::SW);
        bool phase_SE = grid_ptr->get_neighbor_phase(coordinates, DIRECTION::SE);
        return ((phase_N == phase_NE) && (phase_N == phase_NW) && (phase_S == phase_SW) && (phase_S == phase_SE) && (phase_N != phase_S));
    }
}

SUBLATTICE EventRateCalculator::get_sublattice_of_cell(const Coordinates& coordinates) const
{
    int x = coordinates.first;
    int y = coordinates.second;
    bool phase = grid_ptr->get_cell_phase(x, y);
    bool even;
    if (boundary_minus)
    {
        even = (modulo(y, 2) == 0);
    }
    else
    {
        even = (modulo(x + y, 2) == 0);
    }
    if (even == !phase)
    {
        return SUBLATTICE::A;
    }
    else
    {
        return SUBLATTICE::B;
    }
}

double EventRateCalculator::calculate_repulsion_energy_change(const Event& event) const
{
    double repulsion_energy = (boundary_minus ? MINUS_REPULSION : PLUS_REPULSION);
    double energy_change = 0.0;
    for (const auto& coordinates : event)
    {
        bool phase = grid_ptr->get_cell_phase(coordinates);
        bool phase_below = (boundary_minus ? grid_ptr->get_neighbor_phase(coordinates, DIRECTION::SW)
                                           : grid_ptr->get_neighbor_phase(coordinates, DIRECTION::S));
        bool boundary_below = (phase != phase_below);
        bool phase_above_2 = grid_ptr->get_cell_phase(coordinates.first, coordinates.second + 2);
        bool phase_below_2 = grid_ptr->get_cell_phase(coordinates.first, coordinates.second - 2);
        bool phase_closer_to_boundary = (boundary_below ? phase_below_2 : phase_above_2);
        bool phase_further_from_boundary = (boundary_below ? phase_above_2 : phase_below_2);

        if (phase_further_from_boundary != phase)
        {
            energy_change += repulsion_energy;
        }
        if (phase_closer_to_boundary == phase)
        {
            energy_change -= repulsion_energy;
        }
    }
    return energy_change;
}
