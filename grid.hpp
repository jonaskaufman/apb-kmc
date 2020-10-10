#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <utility>
#include <vector>

/// Non-negative modulo
inline int modulo(int a, int b) { return ((a % b) + b) % b; }

// Cardinal directions
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

using PixelGrid = std::vector<std::vector<double>>;

/// Periodic grid of simulation cells with phase values 0 or 1
/// TODO (and sublattice values)
/// Can be staggered rectangular grid (like bricks) or square grid
class SimulationCellGrid
{
public:
    SimulationCellGrid() = delete;

    /// Construct grid of a given width with regions of alternating phase
    //  of heights given by spacings.
    //  An even number of spacings must be given to ensure periodicity
    SimulationCellGrid(const int& width, const std::vector<int>& spacings, bool stagger);

    /// Grid dimensions
    const int width;
    const int height;

    /// Get phase of cell located at x,y
    bool get_cell_phase(const int& x, const int& y) const;

    /// Get indices of a neighbor of cell at x,y specified by direction 
    std::pair<int, int> get_neighbor_indices(const int& x, const int& y, const DIRECTION& direction) const;

    /// Get phase of a neighbor of cell at x,y specified by direction
    bool get_neighbor_phase(const int& x, const int& y, const DIRECTION& direction) const;

    /// Flip phase of cell located at x,y from 0 to 1, or vice versa
    void flip_cell_phase(const int& x, const int& y);

    // TODO: Use typdef or struct for pixel grids?
    /// Get pixel representation of phase grid
    //  If stagger is false, pixel grid is a normal grid of square cells
    //  If stagger is true, cells are doubled along the x direction and staggered along y
    PixelGrid get_phase_pixel_grid() const;

    /// Get pixel representation of grid where the value of each pixel
    //  is equal to the combined vertical distance to the nearest two boundaries
//    std::vector<std::vector<int>> get_spacings_pixel_grid() const;

private:
    /// Whether the grid is staggered
    bool stagger;

    /// Phase values of grid
    std::vector<std::vector<bool>> phase_grid;

    /// Set phase of cell located at x,y
    void set_cell_phase(const int& x, const int& y, const bool& phase);
};

// TODO: Come up with a better name
/// Average pixel values across each row to obtain a vertical profile
std::vector<double> average_horizontal_pixels(const PixelGrid& pixel_grid);

/// TODO is this needed?
/// Print pixel grid to output stream
void print_pixel_grid(const PixelGrid& pixel_grid, std::ostream& stream);

#endif