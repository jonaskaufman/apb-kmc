#include "calculator.hpp"
#include "grid.hpp"
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
    // Cell must be at boundary
    int x = event[0].first;
    int y = event[0].second;
    bool phase_NW = grid.get_neighbor_phase(x, y, DIRECTION::NW);
    bool phase_SW = grid.get_neighbor_phase(x, y, DIRECTION::SW);
    if (phase_NW == phase_SW)
    {
        return false;
    }
    // Boundary must be flat
    bool phase_NE = grid.get_neighbor_phase(x, y, DIRECTION::NE);
    bool phase_SE = grid.get_neighbor_phase(x, y, DIRECTION::SE);
    if ((phase_NW != phase_NE) || (phase_SW != phase_SE))
    {
        return false;
    }
    return true;
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
        // Cell must be at boundary
        int x = coordinates.first;
        int y = coordinates.second;
        bool phase_N = grid.get_neighbor_phase(x, y, DIRECTION::N);
        bool phase_S = grid.get_neighbor_phase(x, y, DIRECTION::S);
        if (phase_N == phase_S)
        {
            return false;
        }
        // Boundary must be flat
        bool phase_NW = grid.get_neighbor_phase(x, y, DIRECTION::NW);
        bool phase_SW = grid.get_neighbor_phase(x, y, DIRECTION::SW);
        bool phase_NE = grid.get_neighbor_phase(x, y, DIRECTION::NE);
        bool phase_SE = grid.get_neighbor_phase(x, y, DIRECTION::SE);
        if ((phase_NW != phase_N) || (phase_NE != phase_N) || (phase_SW != phase_S) || (phase_SE != phase_S))
        {
            return false;
        }
    }
    if (event.size() == 1)
    {
        int x = event[0].first;
        int y = event[0].second;
        bool phase = grid.get_cell_phase(x, y);
        bool phase_W = grid.get_neighbor_phase(x, y, DIRECTION::W);
        bool phase_E = grid.get_neighbor_phase(x, y, DIRECTION::E);
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

//// TODO break up into smaller functions
// double EventRateCalculator::calculate_total_repulsion_energy_change(const Event& event) const
//{
//    double energy_change = 0;
//    for (auto& indices : event)
//    {
//        if (boundary_type == BOUNDARY_TYPE::MINUS)
//        {
//            int x = indices.first;
//            int y = indices.second;
//            energy_change += calculate_repulsion_energy_change(x, y);
//        }
//    }
//    return energy_change;
//}
//
//// TODO: Change to nearest-neighbor boundary only
// double EventRateCalculator::calculate_repulsion_energy_change(int x, int y) const
//{
//    double energy_change = 0;
//    int phase = grid.get_cell_phase(x, y);
//    bool boundary_below = phase != grid.get_cell_phase(x, y - 1);
//    // TODO assert site is at boundary
//
//    bool found_boundary_same = false;
//    bool found_boundary_opposite = false;
//    int d_same = 100;
//    int d_opposite = 100;
//
//    int max_dy = 5;
//    for (int i = 0; i < 2; i++)
//    {
//        for (int dy = 1; dy <= max_dy; dy++)
//        {
//            int x_look = (i == 0 ? (x - 1 * (y % 2) * (dy % 2)) : (x + 1 * ((y + 1) % 2) * (dy % 2)));
//
//            if (found_boundary_same && found_boundary_opposite)
//            {
//                break;
//            }
//            if (!found_boundary_same)
//            {
//                if ((boundary_below && (phase != grid.get_cell_phase(x_look, y + dy))) ||
//                    (!boundary_below && (phase != grid.get_cell_phase(x_look, y - dy))))
//                {
//                    found_boundary_same = true;
//                    d_same = dy;
//                }
//            }
//            if (!found_boundary_opposite)
//            {
//                if ((boundary_below && (phase == grid.get_cell_phase(x_look, y - dy))) ||
//                    (!boundary_below && (phase == grid.get_cell_phase(x_look, y + dy))))
//                {
//                    found_boundary_opposite = true;
//                    d_opposite = dy - 1;
//                }
//            }
//        }
//        double current_energy = 0.5 * (zeta_minus_boundary_energy(d_same) + zeta_minus_boundary_energy(d_opposite));
//        double new_energy = 0.5 * (zeta_minus_boundary_energy(d_same - 1) + zeta_minus_boundary_energy(d_opposite +
//        1)); energy_change += new_energy - current_energy;
//        //            std::cout << d_same << ", " << d_opposite << ", " << new_energy - current_energy <<
//        //            std::endl;
//    }
//    return energy_change;
//}

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

