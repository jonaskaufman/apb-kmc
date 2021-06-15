#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <utility>
#include <vector>

// TODO: Move other aliases here?

// Kinetic event represented by indices of grid cells to flip
using Event = std::vector<std::pair<int, int>>;

using ID = int;

/// Antiphase boundary types (zeta minus or zeta plus)
enum class BOUNDARY_TYPE
{
    MINUS,
    PLUS
};

/// Sublattices on honeycomb network
enum class SUBLATTICE
{
    A,
    B
};

#endif
