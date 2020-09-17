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



void SimulationCellGrid::print_pixel_grid(std::ostream& stream) const
{
    // TODO add staggering properly
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            stream << phase_grid[x][y] << " ";
        }
        stream << std::endl;
    }
}

    int SimulationCellGrid::get_cell_phase(const int& x, const int& y) const
    {
        return phase_grid[x % width][y % height];
    }

    void SimulationCellGrid::flip_cell_phase(const int& x, const int& y)
    {
        int current_phase = get_cell_phase(x, y);
        set_cell_phase(x, y, (current_phase + 1) % 2);
    }

    void SimulationCellGrid::set_cell_phase(const int& x, const int& y, const int& phase)
    {
        phase_grid[x % width][y % height] = phase;
    }

