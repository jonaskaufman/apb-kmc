#include "grid.hpp"
#include <numeric>

SimulationCellGrid::SimulationCellGrid(const int& width, const std::vector<int>& spacings, bool stagger)
    : width{width}, height{std::accumulate(spacings.begin(), spacings.end(), 0)}, stagger{stagger}
{
    phase_grid = std::vector<std::vector<bool>>(width, std::vector<bool>(height));
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

std::vector<std::vector<int>> SimulationCellGrid::get_phase_pixel_grid() const
{
    int pixel_grid_width = (stagger ? 2 * width : width);
    std::vector<std::vector<int>> pixel_grid(pixel_grid_width, std::vector<int>(height));
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            int phase = get_cell_phase(x, y);
            if (stagger)
            {
                pixel_grid[2 * x][y] = phase;
                pixel_grid[modulo(2 * x + 1 - 2 * (y % 2), pixel_grid_width)][y] = phase;
            }
            else
            {
                pixel_grid[x][y] = phase;
            }
        }
    }
    return pixel_grid;
}

std::vector<std::vector<int>> SimulationCellGrid::get_spacings_pixel_grid() const
{
    std::vector<std::vector<int>> phase_pixel_grid = get_phase_pixel_grid();
    int pixel_grid_width = phase_pixel_grid.size();
    int pixel_grid_height = phase_pixel_grid[0].size();
    std::vector<std::vector<int>> spacings_pixel_grid(pixel_grid_width, std::vector<int>(height));
    for (int y = 0; y < pixel_grid_height; y++)
    {
        for (int x = 0; x < pixel_grid_width; x++)
        {
            int phase = phase_pixel_grid[x][y];
            int k_up = -1;
            for (int dy = 0; dy < pixel_grid_height; dy++)
            {
                if (phase != phase_pixel_grid[x][modulo(y + dy, pixel_grid_height)])
                {
                    k_up = dy;
                    break;
                }
            }
            int k_down = -1;
            for (int dy = 0; dy < pixel_grid_height; dy++)
            {
                if (phase != phase_pixel_grid[x][modulo(y - dy, pixel_grid_height)])
                {
                    k_down = dy;
                    break;
                }
            }

            int k = k_up + k_down - 1;
            spacings_pixel_grid[x][y] = k;
        }
    }
    return spacings_pixel_grid;
}

bool SimulationCellGrid::get_cell_phase(const int& x, const int& y) const
{
    return phase_grid[modulo(x, width)][modulo(y, height)];
}

void SimulationCellGrid::flip_cell_phase(const int& x, const int& y)
{
    int current_phase = get_cell_phase(x, y);
    set_cell_phase(x, y, !current_phase);
}

void SimulationCellGrid::set_cell_phase(const int& x, const int& y, const bool& phase)
{
    phase_grid[modulo(x, width)][modulo(y, height)] = phase;
}

std::vector<double> get_horizontal_pixel_averages(const std::vector<std::vector<int>>& pixel_grid)
{
    int pixel_grid_width = pixel_grid.size();
    int pixel_grid_height = pixel_grid[0].size();
    std::vector<double> averages(pixel_grid_height, 0);
    for (int y = 0; y < pixel_grid_height; y++)
    {
        int sum = 0;
        for (int x = 0; x < pixel_grid_width; x++)
        {
            sum += pixel_grid[x][y];
        }
        averages[y] = static_cast<double>(sum) / static_cast<double>(pixel_grid_width);
    }
    return averages;
}

void print_pixel_grid(const std::vector<std::vector<int>>& pixel_grid, std::ostream& stream)
{
    for (int y = 0; y < pixel_grid[0].size(); y++)
    {
        for (int x = 0; x < pixel_grid.size(); x++)
        {
            stream << pixel_grid[x][y] << " ";
        }
        stream << std::endl;
    }
}

