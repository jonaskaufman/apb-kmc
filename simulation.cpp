#include "simulation.hpp"
#include "parameters.hpp"
#include <random>
#include <stdexcept>

/// Boltzmann factor
inline double boltzmann_factor(double barrier, double temperature) { return std::exp(-barrier / (KB * temperature)); }

RandomGenerator::RandomGenerator()
{
    std::random_device device;
    generator = std::mt19937_64(device());
}

int RandomGenerator::sample_integer_range(int maximum_value)
{
    return std::uniform_int_distribution<int>(0, maximum_value)(generator);
}

double RandomGenerator::sample_unit_interval() { return std::uniform_real_distribution<double>(0.0, 1.0)(generator); }

Simulation::Simulation(BOUNDARY_TYPE boundary_type, const SimulationCellGrid& initial_grid, int temperature)
    : boundary_type{boundary_type}, grid{initial_grid}, time{0}, temperature{temperature}
{
    // Check if grid is compatible with boundary type
    bool correct_staggered = (boundary_type == BOUNDARY_TYPE::MINUS);
    if (grid.staggered != correct_staggered)
    {
        throw std::runtime_error("Boundary type not compatible with grid.");
    }
    populate_event_list();
}

void Simulation::step()
{
    Event candidate_event = event_list[random_generator.sample_integer_range(event_list.size() - 1)];
    double rate = calculate_rate(candidate_event);
    if (rate < 0 || rate > 1)
    {
        throw std::runtime_error("Rate is invalid.");
    }
    if (rate > random_generator.sample_unit_interval())
    {
        for (auto& indices : candidate_event)
        {
            grid.flip_cell_phase(indices.first, indices.second);
        }
    }
    double time_step = std::log(1.0 / random_generator.sample_unit_interval()) / event_list.size();
    if (time_step < 0)
    {
        throw std::runtime_error("Time step is invalid.");
    }
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
    // Find the y coordinate of each boundary at x = 0
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

    // zeta plus currently wrong!
    int n_boundaries = boundary_y_origins.size();
    bool minus = (boundary_type == BOUNDARY_TYPE::MINUS);
    int composition_grid_height =
        (minus ? (4 * phase_grid_height + 2 * n_boundaries) : (4 * phase_grid_height - n_boundaries));
    PixelGrid composition_grid(width, std::vector<double>(composition_grid_height, 0.5));
    // Follow each boundary along the x direction and update composition accordingly
    for (int b = 0; b < n_boundaries; b++)
    {
        int y_phase = boundary_y_origins[b];
        int y_comp = (minus ? (4 * y_phase + 2 * b) : (4 * y_phase - b));
        for (int x = 0; x < width; x++)
        {
            if (minus)
            {
                composition_grid[x][modulo(y_comp, composition_grid_height)] = 0.0;
                composition_grid[x][modulo(y_comp + 1, composition_grid_height)] = 0.0;
            }
            else
            {
                composition_grid[x][modulo(y_comp, composition_grid_height)] = 1.0;
            }
            int x_next = modulo(x + 1, width);
            int y_phase_previous = modulo(y_phase - 1, phase_grid_height);
            int y_phase_current = modulo(y_phase, phase_grid_height);
            if (phase_grid[x][y_phase_current] != phase_grid[x_next][y_phase_current])
            {
                y_phase += 1;
                y_comp += 4;
            }
            else if (phase_grid[x][y_phase_previous] != phase_grid[x_next][y_phase_previous])
            {
                y_phase -= 1;
                y_comp -= 4;
            }
        }
    }
    return composition_grid;
}

std::vector<double> Simulation::get_average_composition_profile() const
{
    return average_pixels_horizontally(get_composition_pixel_grid());
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
double Simulation::calculate_repulsion_energy_change(int x, int y) const
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

SUBLATTICE Simulation::get_sublattice_of_cell(int x, int y)
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

