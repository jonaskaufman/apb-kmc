#include "calculator.hpp"
#include <cmath>
#include <stdexcept>

inline double boltzmann_factor(double barrier, double temperature) { return std::exp(-barrier / (KB * temperature)); }

EventRateCalculator::EventRateCalculator(BOUNDARY_TYPE boundary_type, double temperature)
    : boundary_type{boundary_type}, temperature{temperature}
{
}

double EventRateCalculator::calculate_rate(const Event& event, const SimulationCellGrid& grid) const
{
    if (!is_valid_event(event, grid))
    {
        return 0.0;
    }
    else
    {
        double barrier = calculate_barrier(event, grid);
        return boltzmann_factor(barrier, temperature);
    }
}

bool EventRateCalculator::is_valid_event(const Event& event, const SimulationCellGrid& grid) const
{
    if (boundary_type == BOUNDARY_TYPE::MINUS)
    {
        return is_valid_event_zeta_minus(event, grid);
    }
    else
    {
        return is_valid_event_zeta_plus(event, grid);
    }
}

bool EventRateCalculator::is_valid_event_zeta_minus(const Event& event, const SimulationCellGrid& grid) const
{
    // Only single cell events are allowed
    if (event.size() != 1)
    {
        return false;
    }
    // Cell must be at a flat boundary and not be the only cell separating two boundaries
    return at_valid_boundary(event[0], grid);
}

bool EventRateCalculator::is_valid_event_zeta_plus(const Event& event, const SimulationCellGrid& grid) const
{
    // Only single or double cell events are allowed
    if ((event.size() != 1) && (event.size() != 2))
    {
        return false;
    }
    // Cell(s) must be at a flat boundary and not be the only cell separating two boundaries
    for (const auto& coordinates : event)
    {
        if (!at_valid_boundary(coordinates, grid))
        {
            return false;
        }
    }
    // Additional checks, mostly ensuring no single cell kinks are introduced
    if (!passes_additional_checks_zeta_plus(event, grid))
    {
        return false;
    }
    return true;
}

bool EventRateCalculator::passes_additional_checks_zeta_plus(const Event& event, const SimulationCellGrid& grid) const
{
    if (event.size() == 1)
    {
        bool phase = grid.get_cell_phase(event[0]);
        bool phase_W = grid.get_neighbor_phase(event[0], DIRECTION::W);
        bool phase_E = grid.get_neighbor_phase(event[0], DIRECTION::E);
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
            bool phase_X_X = grid.get_neighbor_phase(grid.get_neighbor_indices(event[0], direction_X), direction_X);
            if (phase_X_X != phase_X)
            {
                return false;
            }
        }
    }
    else if (event.size() == 2)
    {
        // Event cells must be in expected order
        if (grid.get_neighbor_indices(event[0], DIRECTION::E) != event[1])
        {
            throw std::runtime_error("Cells for double atom event are in wrong order or not horizontally adjacent.");
        }
        bool phase_0 = grid.get_cell_phase(event[0]);
        bool phase_1 = grid.get_cell_phase(event[1]);
        // Event cells must have same phase
        if (phase_0 != phase_1)
        {
            return false;
        }
        // Neighbors must have same phase
        bool phase_0_W = grid.get_neighbor_phase(event[0], DIRECTION::W);
        bool phase_1_E = grid.get_neighbor_phase(event[1], DIRECTION::E);
        if (phase_0_W != phase_1_E)
        {
            return false;
        }
        // No single cell kinks should be created
        if (phase_0 == phase_0_W)
        {
            bool phase_0_W_W = grid.get_neighbor_phase(grid.get_neighbor_indices(event[0], DIRECTION::W), DIRECTION::W);
            bool phase_1_E_E = grid.get_neighbor_phase(grid.get_neighbor_indices(event[1], DIRECTION::E), DIRECTION::E);
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

double EventRateCalculator::calculate_barrier(const Event& event, const SimulationCellGrid& grid) const
{
    if (boundary_type == BOUNDARY_TYPE::MINUS)
    {
        return calculate_barrier_zeta_minus(event, grid);
    }
    else
    {
        return calculate_barrier_zeta_plus(event, grid);
    }
}

double EventRateCalculator::calculate_barrier_zeta_minus(const Event& event, const SimulationCellGrid& grid) const
{
    bool phase_W = grid.get_neighbor_phase(event[0], DIRECTION::W);
    bool phase_E = grid.get_neighbor_phase(event[0], DIRECTION::E);
    if (phase_W == phase_E)
    {
        return MINUS_KINK_FORM;
    }
    else
    {
        return MINUS_KINK_MOVE;
    }
}

double EventRateCalculator::calculate_barrier_zeta_plus(const Event& event, const SimulationCellGrid& grid) const
{
    if (event.size() == 2)
    {
        bool phase = grid.get_cell_phase(event[0]);
        bool phase_sides = grid.get_neighbor_phase(event[0], DIRECTION::W);
        SUBLATTICE sublattice_0 = get_sublattice_of_cell(event[0], grid);
        if (phase == phase_sides)
        {
            return ((sublattice_0 == SUBLATTICE::B) ? PLUS_KINK_FORM_I : PLUS_KINK_FORM_II);
        }
        else
        {
            return ((sublattice_0 == SUBLATTICE::A) ? PLUS_KINK_DESTROY_I : PLUS_KINK_DESTROY_II);
        }
    }
    else
    {
        bool phase_S = grid.get_neighbor_phase(event[0], DIRECTION::S);
        bool phase_E = grid.get_neighbor_phase(event[0], DIRECTION::E);
        SUBLATTICE sublattice_E = get_sublattice_of_cell(grid.get_neighbor_indices(event[0], DIRECTION::E), grid);
        if (phase_S == phase_E)
        {
            return ((sublattice_E == SUBLATTICE::A) ? PLUS_KINK_MOVE_I : PLUS_KINK_MOVE_II);
        }
        else
        {
            return ((sublattice_E == SUBLATTICE::A) ? PLUS_KINK_MOVE_III : PLUS_KINK_MOVE_IV);
        }
    }
}

bool EventRateCalculator::at_valid_boundary(const std::pair<int, int>& coordinates,
                                            const SimulationCellGrid& grid) const

{
    if (boundary_type == BOUNDARY_TYPE::MINUS)
    {
        bool phase_NW = grid.get_neighbor_phase(coordinates, DIRECTION::NW);
        bool phase_NE = grid.get_neighbor_phase(coordinates, DIRECTION::NE);
        bool phase_SW = grid.get_neighbor_phase(coordinates, DIRECTION::SW);
        bool phase_SE = grid.get_neighbor_phase(coordinates, DIRECTION::SE);
        return ((phase_NW == phase_NE) && (phase_SW == phase_SE) && (phase_NW != phase_SW));
    }
    else
    {
        bool phase_N = grid.get_neighbor_phase(coordinates, DIRECTION::N);
        bool phase_NW = grid.get_neighbor_phase(coordinates, DIRECTION::NW);
        bool phase_NE = grid.get_neighbor_phase(coordinates, DIRECTION::NE);
        bool phase_S = grid.get_neighbor_phase(coordinates, DIRECTION::S);
        bool phase_SW = grid.get_neighbor_phase(coordinates, DIRECTION::SW);
        bool phase_SE = grid.get_neighbor_phase(coordinates, DIRECTION::SE);
        return ((phase_N == phase_NE) && (phase_N == phase_NW) && (phase_S == phase_SW) && (phase_S == phase_SE) &&
                (phase_N != phase_S));
    }
}

SUBLATTICE EventRateCalculator::get_sublattice_of_cell(const std::pair<int, int>& coordinates,
                                                       const SimulationCellGrid& grid) const
{
    int x = coordinates.first;
    int y = coordinates.second;
    bool phase = grid.get_cell_phase(x, y);
    bool even;
    if (boundary_type == BOUNDARY_TYPE::MINUS)
    {
        even = (modulo(y, 2) == 0);
    }
    else // boundary_type == BOUNDARY_TYPE::PLUS
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

