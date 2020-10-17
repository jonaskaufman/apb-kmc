#include "simulation.hpp"
#include "parameters.hpp"
#include <random>

/// Boltzmann factor
inline double boltzmann_factor(double barrier, double temperature) { return std::exp(-barrier / (KB * temperature)); }

RandomGenerator::RandomGenerator()
{
    std::random_device device;
    generator = std::mt19937_64(device());
}

int RandomGenerator::sample_integer_range(const int& maximum_value)
{
    return std::uniform_int_distribution<int>(0, maximum_value)(generator);
}

double RandomGenerator::sample_unit_interval() { return std::uniform_real_distribution<double>(0.0, 1.0)(generator); }

Simulation::Simulation(const BOUNDARY_TYPE& boundary_type,
                       const SimulationCellGrid& initial_grid,
                       const int& temperature)
    : boundary_type{boundary_type}, grid{initial_grid}, time{0}, temperature{temperature}
{
    // TODO check stagger compatible with boundary type
    populate_event_list();
}

void Simulation::step()
{
    Event candidate_event = event_list[random_generator.sample_integer_range(event_list.size() - 1)];
    double rate = calculate_rate(candidate_event);
    if (rate > random_generator.sample_unit_interval())
    {
        for (auto& indices : candidate_event)
        {
            grid.flip_cell_phase(indices.first, indices.second);
        }
    }
    double time_step = std::log(1.0 / random_generator.sample_unit_interval()) / event_list.size();
    time += time_step;
}

void Simulation::pass()
{
    for (int i = 0; i < event_list.size(); i++)
    {
        step();
    }
}

PixelGrid Simulation::get_phase_pixel_grid() const { return grid.get_phase_pixel_grid(); }

PixelGrid Simulation::get_composition_pixel_grid() const
{
    const PixelGrid phase_grid = get_phase_pixel_grid();
    int width = phase_grid.size();
    int phase_grid_height = phase_grid[0].size();
    std::vector<int> boundary_y_origins;
    for (int y = 0; y < phase_grid_height; y++)
    {
        double current_phase = phase_grid[0][y];
        double previous_phase = phase_grid[0][modulo(y - 1, phase_grid_height)];
        if (current_phase != previous_phase)
        {
            boundary_y_origins.push_back(y);
        }
    }
    int n_boundaries = boundary_y_origins.size();
    // TODO fix for zeta plus boundaries too
    int composition_grid_height = 2 * phase_grid_height + n_boundaries;
    PixelGrid composition_grid(width, std::vector<double>(composition_grid_height, 0.5));
    for (int b = 0; b < n_boundaries; b++)
    {
        int y_p = boundary_y_origins[b];
        int y_c = 2 * y_p + b;
        for (int x = 0; x < width; x++)
        {
            composition_grid[x][modulo(y_c, composition_grid_height)] = 0.0;
            int x_next = modulo(x + 1, width);
            int y_previous = modulo(y_p - 1, phase_grid_height);
            int y_current = y_p;
            if (phase_grid[x][y_current] != phase_grid[x_next][y_current])
            {
                y_p = modulo(y_p + 1, phase_grid_height);
                y_c += 2;
            }
            else if (phase_grid[x][y_previous] != phase_grid[x_next][y_previous])
            {
                y_p = modulo(y_p - 1, phase_grid_height);
                y_c -= 2;
            }
        }
    }
    return composition_grid;
}

std::vector<double> Simulation::average_horizontal_composition_pixels() const
{
    return average_horizontal_pixels(get_composition_pixel_grid());
}

void Simulation::populate_event_list()
{
    event_list.clear();
    for (int y = 0; y < grid.height; y++)
    {
        for (int x = 0; x < grid.width; x++)
        {
            event_list.push_back({std::make_pair(x, y)});
            if (boundary_type == BOUNDARY_TYPE::PLUS)
            {
                event_list.push_back({std::make_pair(x, y), std::make_pair(modulo(x + 1, grid.width), y)});
            }
        }
    }
}

// TODO break up into: Check valid, calculate barrier functions?
double Simulation::calculate_rate(const Event& event) const
{
    double rate = 0;
    if (boundary_type == BOUNDARY_TYPE::MINUS)
    {
        // TODO check that event is single-atom hop
        int x = event[0].first;
        int y = event[0].second;
        int phase = grid.get_cell_phase(x, y);

        // Above/below neighbors are at x, x+1 if y is even and x-1, x if y is odd
        int x_left = x - 1 * (y % 2);
        int x_right = x + 1 * ((y + 1) % 2);
        int up_left_phase = grid.get_cell_phase(x_left, y + 1);
        int up_right_phase = grid.get_cell_phase(x_right, y + 1);
        int down_left_phase = grid.get_cell_phase(x_left, y - 1);
        int down_right_phase = grid.get_cell_phase(x_right, y - 1);

        //        bool flat_check = () && ();

        bool flat_check = (up_left_phase == up_right_phase) && (down_left_phase == down_right_phase);
        bool boundary_check = (phase != up_left_phase) ^ (phase != down_left_phase);
        if (flat_check && boundary_check)
        {
            double barrier = ZETA_MINUS_BARRIER + 0.5 * calculate_total_repulsion_energy_change(
                                                            event); // TODO calculate barrier in separate function?
            rate = boltzmann_factor(barrier, temperature);
        }
    }
    return rate;
}

// TODO break up into smaller functions
double Simulation::calculate_total_repulsion_energy_change(const Event& event) const
{
    double energy_change = 0;
    for (auto& indices : event)
    {
        if (boundary_type == BOUNDARY_TYPE::MINUS)
        {
            int x = indices.first;
            int y = indices.second;
            energy_change += calculate_repulsion_energy_change(x, y);
        }
    }
    return energy_change;
}

// TODO: Change to nearest-neighbor boundary only
double Simulation::calculate_repulsion_energy_change(const int& x, const int& y) const
{
    double energy_change = 0;
    int phase = grid.get_cell_phase(x, y);
    bool boundary_below = phase != grid.get_cell_phase(x, y - 1);
    // TODO assert site is at boundary

    bool found_boundary_same = false;
    bool found_boundary_opposite = false;
    int d_same = 100;
    int d_opposite = 100;

    int max_dy = 5;
    for (int i = 0; i < 2; i++)
    {
        for (int dy = 1; dy <= max_dy; dy++)
        {
            int x_look = (i == 0 ? (x - 1 * (y % 2) * (dy % 2)) : (x + 1 * ((y + 1) % 2) * (dy % 2)));

            if (found_boundary_same && found_boundary_opposite)
            {
                break;
            }
            if (!found_boundary_same)
            {
                if ((boundary_below && (phase != grid.get_cell_phase(x_look, y + dy))) ||
                    (!boundary_below && (phase != grid.get_cell_phase(x_look, y - dy))))
                {
                    found_boundary_same = true;
                    d_same = dy;
                }
            }
            if (!found_boundary_opposite)
            {
                if ((boundary_below && (phase == grid.get_cell_phase(x_look, y - dy))) ||
                    (!boundary_below && (phase == grid.get_cell_phase(x_look, y + dy))))
                {
                    found_boundary_opposite = true;
                    d_opposite = dy - 1;
                }
            }
        }
        double current_energy = 0.5 * (zeta_minus_boundary_energy(d_same) + zeta_minus_boundary_energy(d_opposite));
        double new_energy = 0.5 * (zeta_minus_boundary_energy(d_same - 1) + zeta_minus_boundary_energy(d_opposite + 1));
        energy_change += new_energy - current_energy;
        //            std::cout << d_same << ", " << d_opposite << ", " << new_energy - current_energy <<
        //            std::endl;
    }
    return energy_change;
}

SUBLATTICE Simulation::get_sublattice_of_cell(const int& x, const int& y)
{
    bool phase = grid.get_cell_phase(x, y);
    bool odd;
    if (boundary_type == BOUNDARY_TYPE::MINUS)
    {
        odd = (modulo(y, 2) == 0);
    }
    else // boundary_type == BOUNDARY_TYPE::PLUS
    {
        odd = (modulo(x + y, 2) == 0);
    }
    if (odd == phase)
    {
        return SUBLATTICE::A;
    }
    else
    {
        return SUBLATTICE::B;
    }
}

