#include "simulation.hpp"
#include "parameters.hpp"
#include <iostream>
#include <numeric>
#include <random>

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
    event_index_distribution = std::uniform_int_distribution<int>(0, event_list.size()-1);
    uniform_unit_interval_distribution = std::uniform_real_distribution<double>(0.0, 1.0);
    std::random_device device;
    generator = std::mt19937_64(device());
}


void Simulation::step()
{
    // select random number in size of event list 
    std::cout << event_index_distribution(generator) << std::endl;
    // determine rate of event

    // select number to determine acceptance
    std::cout << uniform_unit_interval_distribution(generator) << std::endl;
    // select number to update time
    double time_step = 0;
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
