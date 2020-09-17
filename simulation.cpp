#include "simulation.hpp"
#include "parameters.hpp"
#include <iostream>
#include <numeric>

SimulationCellGrid::SimulationCellGrid(const int& width, const std::vector<int>& spacings, bool stagger)
    : width{width}, stagger{stagger}
{
    height = std::accumulate(spacings.begin(), spacings.end(), 0);
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

