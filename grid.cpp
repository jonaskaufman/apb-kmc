#include "grid.hpp"
#include <numeric>
#include <stdexcept>

void print_pixel_grid(const PixelGrid& pixel_grid, std::ostream& output_stream)
{
    for (int y = 0; y < pixel_grid[0].size(); y++)
    {
        for (int x = 0; x < pixel_grid.size(); x++)
        {
            output_stream << pixel_grid[x][y] << " ";
        }
        output_stream << std::endl;
    }
}

std::vector<double> average_pixels_horizontally(const PixelGrid& pixel_grid)
{
    int pixel_grid_width = pixel_grid.size();
    int pixel_grid_height = pixel_grid[0].size();
    std::vector<double> averages(pixel_grid_height, 0);
    for (int y = 0; y < pixel_grid_height; y++)
    {
        double total = 0.0;
        for (int x = 0; x < pixel_grid_width; x++)
        {
            total += pixel_grid[x][y];
        }
        averages[y] = total / pixel_grid_width;
    }
    return averages;
}

SimulationCellGrid::SimulationCellGrid(int width, const std::vector<int>& block_heights, bool staggered)
    : width{width}, height{std::accumulate(block_heights.begin(), block_heights.end(), 0)}, staggered{staggered}
{
    // Checks for periodicity
    if (width % 2 != 0 || height % 2 != 0)
    {
        throw std::runtime_error("Grid dimensions must be even to achieve periodicity.");
    }
    if (block_heights.size() % 2 != 0)
    {
        throw std::runtime_error("Number of blocks must be even to achieve periodicity.");
    }
    // Initialize phase grid
    phase_grid = std::vector<std::vector<bool>>(width, std::vector<bool>(height, false));
    int y = 0;
    bool flipped = false;
    for (const auto& spacing : block_heights)
    {
        for (int y_s = 0; y_s < spacing; y_s++)
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

bool SimulationCellGrid::get_cell_phase(int x, int y) const { return phase_grid[modulo(x, width)][modulo(y, height)]; }

bool SimulationCellGrid::get_cell_phase(const Coordinates& coordinates) const
{
    return get_cell_phase(coordinates.first, coordinates.second);
}

Coordinates SimulationCellGrid::get_neighbor_indices(int x, int y, DIRECTION neighbor_direction) const
{
    // Check that direction is valid
    if (staggered && (neighbor_direction == DIRECTION::N || neighbor_direction == DIRECTION::S))
    {
        throw std::runtime_error("N and S are invalid neighbor directions for a staggered grid.");
    }
    // Adjust x coordinate, if direction has E/W component
    int dx = 0;
    if (neighbor_direction == DIRECTION::E || neighbor_direction == DIRECTION::NE ||
        neighbor_direction == DIRECTION::SE)
    {
        if (staggered && (neighbor_direction != DIRECTION::E))
        {
            dx = 1 * ((y + 1) % 2);
        }
        else
        {
            dx = 1;
        }
    }
    else if (neighbor_direction == DIRECTION::W || neighbor_direction == DIRECTION::NW ||
             neighbor_direction == DIRECTION::SW)
    {
        if (staggered && (neighbor_direction != DIRECTION::W))
        {
            dx = -1 * (y % 2);
        }
        else
        {
            dx = -1;
        }
    }
    // Adjust y coordinate if direction has N/S component
    int dy = 0;
    if (neighbor_direction == DIRECTION::N || neighbor_direction == DIRECTION::NE ||
        neighbor_direction == DIRECTION::NW)
    {
        dy = 1;
    }
    else if (neighbor_direction == DIRECTION::S || neighbor_direction == DIRECTION::SE ||
             neighbor_direction == DIRECTION::SW)
    {
        dy = -1;
    }
    // Check that indices make sense
    if ((dx == 0 && dy == 0) || dx < -1 || dx > 1 || dy < -1 || dy > 1)
    {
        throw std::runtime_error("Neighbor indices do not correspond to a neighbor.");
    }
    return std::make_pair(modulo(x + dx, width), modulo(y + dy, height));
}

Coordinates SimulationCellGrid::get_neighbor_indices(const Coordinates& coordinates,
                                                             DIRECTION neighbor_direction) const
{
    return get_neighbor_indices(coordinates.first, coordinates.second, neighbor_direction);
}

bool SimulationCellGrid::get_neighbor_phase(int x, int y, DIRECTION neighbor_direction) const
{
    Coordinates neighbor_indices = get_neighbor_indices(x, y, neighbor_direction);
    return get_cell_phase(neighbor_indices.first, neighbor_indices.second);
}

bool SimulationCellGrid::get_neighbor_phase(const Coordinates& coordinates, DIRECTION neighbor_direction) const
{
    return get_neighbor_phase(coordinates.first, coordinates.second, neighbor_direction);
}

void SimulationCellGrid::flip_cell_phase(int x, int y)
{
    int current_phase = get_cell_phase(x, y);
    set_cell_phase(x, y, !current_phase);
}

void SimulationCellGrid::set_cell_phase(int x, int y, bool phase)
{
    phase_grid[modulo(x, width)][modulo(y, height)] = phase;
}

PixelGrid SimulationCellGrid::get_phase_pixel_grid() const
{
    // Double width always
    int pixel_grid_width = 2 * width;
    PixelGrid pixel_grid(pixel_grid_width, std::vector<double>(height, -1.0));
    for (int y = 0; y < height; y++)
    {
        // Set values for current row
        for (int x = 0; x < width; x++)
        {
            // Each cell becomes two pixels
            int phase = get_cell_phase(x, y);
            int x_next = staggered ? (2 * x + 1 - 2 * (y % 2)) : (2 * x + 1);
            pixel_grid[2 * x][y] = phase;
            pixel_grid[modulo(x_next, pixel_grid_width)][y] = phase;
        }
        // Check that all values have been set
        for (int x = 0; x < pixel_grid_width; x++)
        {
            if (pixel_grid[x][y] < 0)
            {
                throw std::runtime_error("Some pixel grid values were not set properly");
            }
        }
    }
    return pixel_grid;
}
