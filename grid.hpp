#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <utility>
#include <vector>

/// Calculates the non-negative value of a modulo b
inline int modulo(int a, int b) { return ((a % b) + b) % b; }

/// Cardinal directions
enum class DIRECTION
{
    N,
    S,
    E,
    W,
    NE,
    NW,
    SE,
    SW
};

/// Two-dimensional square grid of pixel values (x is first dimension, y is second dimension)
using PixelGrid = std::vector<std::vector<double>>;

/// Prints a pixel grid to an output stream
void print_pixel_grid(const PixelGrid& pixel_grid, std::ostream& output_stream);

/// Averages each row of a pixel grid and returns the resulting profile
std::vector<double> average_pixels_horizontally(const PixelGrid& pixel_grid);

/**
 * Periodic two-dimensional grid of cells with phase values false (0) or true (1)
 *
 * Grid can either be a square grid (like graph paper), or a staggered grid
 * in which every other row is offset in the x direction by half a cell width (like a brick wall).
 */
class SimulationCellGrid
{
public:
    SimulationCellGrid() = delete;

    /// Constructs a grid with blocks of alternating phase in the y direction
    SimulationCellGrid(int width, const std::vector<int>& block_heights, bool staggered);

    /// Grid dimensions
    const int width;
    const int height;

    /// Whether the grid is staggered
    const bool staggered;

    /// Returns the phase of the cell at (x, y)
    bool get_cell_phase(int x, int y) const;

    /// Returns the coordinates of the given neighbor of the cell at (x, y)
    std::pair<int, int> get_neighbor_indices(int x, int y, DIRECTION neighbor_direction) const;

    /// Returns the phase of the given neighbor of the cell at (x, y)
    bool get_neighbor_phase(int x, int y, DIRECTION neighbor_direction) const;

    /// Flips the phase of the cell at (x, y) to the opposite of its current value
    void flip_cell_phase(int x, int y);

    /// Returns a pixel representation of the grid with phases 0 or 1. Width is always doubled.
    PixelGrid get_phase_pixel_grid() const;

private:
    /// Phase values of the grid
    std::vector<std::vector<bool>> phase_grid;

    /// Sets the phase value of the cell at (x, y)
    void set_cell_phase(int x, int y, bool phase);
};

#endif
