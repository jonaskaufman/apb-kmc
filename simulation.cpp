#include "simulation.hpp"
#include "parameters.hpp"
#include <iostream>
#include <numeric>
#include <random>

/// Positive modulo
inline int modulo(int a, int b) { return ((a % b) + b) % b; }

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

SimulationCellGrid::SimulationCellGrid(const int& width, const std::vector<int>& spacings, bool stagger)
    : width{width}, height{std::accumulate(spacings.begin(), spacings.end(), 0)}, stagger{stagger}
{
    phase_grid = std::vector<std::vector<int>>(width, std::vector<int>(height));
    int y = 0;
    bool flipped = false;
    for (auto& spacing : spacings)
    {
        for (int s = 0; s < spacing; s++)
        {
            for (int x = 0; x < width; x++)
            {
                if (flipped)
                {
                    flip_cell_phase(x, y);
                }
            }
            y++;
        }
        flipped = !flipped;
    }
}

std::vector<std::vector<int>> SimulationCellGrid::get_pixel_grid() const
{
    if (stagger)
    {
        int pixel_grid_width = 2 * width;
        std::vector<std::vector<int>> pixel_grid(pixel_grid_width, std::vector<int>(height));
        for (int y = 0; y < height; y++)
        {
            for (int x = 0; x < width; x++)
            {
                int phase = get_cell_phase(x, y);
                pixel_grid[2 * x][y] = phase;
                pixel_grid[modulo(2 * x + 1 - 2 * (y % 2), pixel_grid_width)][y] = phase;
            }
        }
        return pixel_grid;
    }
    else
    {
        return phase_grid;
    }
}

void SimulationCellGrid::print_pixel_grid(std::ostream& stream) const
{
    std::vector<std::vector<int>> pixel_grid = get_pixel_grid();
    for (int y = 0; y < pixel_grid[0].size(); y++)
    {
        for (int x = 0; x < pixel_grid.size(); x++)
        {
            stream << pixel_grid[x][y] << " ";
        }
        stream << std::endl;
    }
}

int SimulationCellGrid::get_cell_phase(const int& x, const int& y) const
{
    return phase_grid[modulo(x, width)][modulo(y, height)];
}

void SimulationCellGrid::flip_cell_phase(const int& x, const int& y)
{
    int current_phase = get_cell_phase(x, y);
    set_cell_phase(x, y, (current_phase + 1) % 2);
}

void SimulationCellGrid::set_cell_phase(const int& x, const int& y, const int& phase)
{
    phase_grid[modulo(x, width)][modulo(y, height)] = phase;
}

Simulation::Simulation(const BOUNDARY_TYPE& boundary_type,
                       const SimulationCellGrid& initial_grid,
                       const int& temperature)
    : boundary_type{boundary_type}, grid{initial_grid}, time{0}, temperature{temperature}
{
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
    double time_step =
        std::log(1.0 / random_generator.sample_unit_interval()) / event_list.size(); // TODO check for int division here
    time += time_step;
}

void Simulation::print_grid(std::ostream& stream) const { grid.print_pixel_grid(stream); }

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

double Simulation::calculate_rate(const Event& event) const
{
    double rate = 0;
    if (boundary_type == BOUNDARY_TYPE::MINUS)
    {
        // TODO check that event is single-atom hop
        int x = event[0].first;
        int y = event[0].second;
        int phase = grid.get_cell_phase(x, y);
        int above_phase_0 = grid.get_cell_phase(x, y + 1);
        int above_phase_1 = grid.get_cell_phase(x + 1 - 2 * (y % 2), y + 1);
        int below_phase_0 = grid.get_cell_phase(x, y - 1);
        int below_phase_1 = grid.get_cell_phase(x + 1 - 2 * (y % 2), y - 1);
        if (((above_phase_0 == above_phase_1) && (below_phase_0 == below_phase_1)) &&
            ((phase != above_phase_0) ^ (phase != below_phase_0)))
        {
            rate = boltzmann_factor(ZETA_MINUS_BARRIER, temperature);
        }
    }
    return rate;
}
