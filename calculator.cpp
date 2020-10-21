#include "calculator.hpp"
#include <cmath>

inline double boltzmann_factor(double barrier, double temperature) { return std::exp(-barrier / (KB * temperature)); }

EventRateCalculator::EventRateCalculator(BOUNDARY_TYPE boundary_type, double temperature)
    : boundary_type{boundary_type}, temperature{temperature}
{
}

double EventRateCalculator::calculate_rate(const Event& event, const SimulationCellGrid& grid) const
{
    if (!is_valid_event(event, grid))
    {
        return 0;
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
    // Only single atom/cell events are allowed
    if (event.size() != 1)
    {
        return false;
    }
    // Cell must be at a flat boundary, and not be the only cell separating two boundaries
    bool above = boundary_above(event[0], grid);
    bool below = boundary_below(event[0], grid);
    return (above ^ below);
}

bool EventRateCalculator::is_valid_event_zeta_plus(const Event& event, const SimulationCellGrid& grid) const
{
    // Only single or double atom/cell events are allowed
    if ((event.size() != 1) && (event.size() != 2))
    {
        return false;
    }
    for (const auto& coordinates : event)
    {
        // Cell must be at a flat boundary, and not be the only cell separating two boundaries
        bool above = boundary_above(event[0], grid);
        bool below = boundary_below(event[0], grid);
        if (!(above ^ below))
        {
            return false;
        }
    }
    // TODO
    if (event.size() == 1)
    {
        bool phase = grid.get_cell_phase(event[0]);
        bool phase_W = grid.get_neighbor_phase(event[0], DIRECTION::W);
        bool phase_E = grid.get_neighbor_phase(event[0], DIRECTION::E);
        if (phase_W == phase_E)
        {
            if (phase != phase_W)
            {
                throw std::runtime_error("Single atom kink detected for zeta plus, which is not allowed.");
            }
            return false;
        }
        // TODO check about single kink formation
    }
    else
    {
        // Cells must have the same phase
        int x_0 = event[0].first;
        int x_1 = event[1].first;
        int y_0 = event[0].second;
        int y_1 = event[1].second;
        if ((y_0 != y_1) || (x_1 != modulo(x_0 + 1, grid.width)))
        {
            throw std::runtime_error("Cells for double atom event are not adjacent, or in wrong order.");
        }
        int y = y_0;
        bool phase_0 = grid.get_cell_phase(x_0, y);
        bool phase_1 = grid.get_cell_phase(x_1, y);
        if (phase_0 != phase_1)
        {
            return false;
        }
        bool phase_W_0 = grid.get_neighbor_phase(x_0, y, DIRECTION::W);
        bool phase_E_1 = grid.get_neighbor_phase(x_1, y, DIRECTION::E);
        if (phase_W_0 != phase_E_1)
        {
            return false;
        }
        bool phase_W_W_0 = grid.get_neighbor_phase(grid.get_neighbor_indices(x_0, y, DIRECTION::W), DIRECTION::W);
        bool phase_E_E_1 = grid.get_neighbor_phase(grid.get_neighbor_indices(x_1, y, DIRECTION::E), DIRECTION::E);
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

bool EventRateCalculator::boundary_above(const std::pair<int, int>& coordinates, const SimulationCellGrid& grid) const
{
    bool phase = grid.get_cell_phase(coordinates);
    if (boundary_type == BOUNDARY_TYPE::MINUS)
    {
        bool phase_NW = grid.get_neighbor_phase(coordinates, DIRECTION::NW);
        bool phase_NE = grid.get_neighbor_phase(coordinates, DIRECTION::NE);
        return ((phase_NW == phase_NE) && (phase != phase_NW));
    }
    else
    {
        bool phase_N = grid.get_neighbor_phase(coordinates, DIRECTION::N);
        bool phase_NW = grid.get_neighbor_phase(coordinates, DIRECTION::NW);
        bool phase_NE = grid.get_neighbor_phase(coordinates, DIRECTION::NE);
        return ((phase_N == phase_NW) && (phase_N == phase_NE) && (phase != phase_N));
    }
}

bool EventRateCalculator::boundary_below(const std::pair<int, int>& coordinates, const SimulationCellGrid& grid) const
{
    bool phase = grid.get_cell_phase(coordinates);
    if (boundary_type == BOUNDARY_TYPE::MINUS)
    {
        bool phase_SW = grid.get_neighbor_phase(coordinates, DIRECTION::SW);
        bool phase_SE = grid.get_neighbor_phase(coordinates, DIRECTION::SE);
        return ((phase_SW == phase_SE) && (phase != phase_SW));
    }
    else
    {
        bool phase_S = grid.get_neighbor_phase(coordinates, DIRECTION::S);
        bool phase_SW = grid.get_neighbor_phase(coordinates, DIRECTION::SW);
        bool phase_SE = grid.get_neighbor_phase(coordinates, DIRECTION::SE);
        return ((phase_S == phase_SW) && (phase_S == phase_SE) && (phase != phase_S));
    }
}

SUBLATTICE EventRateCalculator::get_sublattice_of_cell(int x, int y, const SimulationCellGrid& grid)
{
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

