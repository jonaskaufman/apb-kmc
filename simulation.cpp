#include "simulation.hpp"
#include <random>
#include <stdexcept>

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

Simulation::Simulation(BOUNDARY_TYPE boundary_type, const SimulationCellGrid& initial_grid, double temperature)
    : boundary_type{boundary_type},
      grid{initial_grid},
      time{0},
      temperature{temperature},
      rate_calculator{boundary_type, temperature}
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
    if (rate != 0.0)
    {
        if (rate > random_generator.sample_unit_interval())
        {
            for (auto& indices : candidate_event)
            {
                grid.flip_cell_phase(indices.first, indices.second);
            }
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

PixelGrid Simulation::get_phase_pixel_grid() const
{
    const PixelGrid phase_grid = grid.get_phase_pixel_grid();
    int width = phase_grid.size();
    int height = phase_grid[0].size();
    int scaled_height = 4 * height;
    PixelGrid scaled_phase_grid(width, std::vector<double>(scaled_height, 0));
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            for (int dy = 0; dy < 4; dy++)
            {
                scaled_phase_grid[x][4 * y + dy] = phase_grid[x][y];
            }
        }
    }
    return scaled_phase_grid;
}

PixelGrid Simulation::get_composition_pixel_grid() const
{
    const PixelGrid phase_grid = grid.get_phase_pixel_grid();
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

double Simulation::calculate_rate(const Event& event) const { return rate_calculator.calculate_rate(event, grid); }
