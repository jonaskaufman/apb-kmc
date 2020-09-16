#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>
#include <utility>
#include <vector>

/// Periodic grid of simulation cells
class SimulationCellGrid
{
public:
    SimulationCellGrid() = delete;
    SimulationCellGrid(const int& width, const std::vector<int>& spacings);
    int get_cell_phase(const int& x, const int& y) const;
    void flip_cell_phase(const int& x, const int& y);
    const std::vector<std::vector<int>>& get_pixel_grid() const; // reference?
    void print_grid(std::ostream& stream) const;

private:
    int width;
    int height;
    bool stagger;
    std::vector<std::vector<int>> phase_grid;
    void set_cell_phase(int x, int y, int phase);
};

/// Kinetic events are vectors of indices representing sites to flip 
using Event = std::vector<std::pair<int, int>>;

/// Monte Carlo simulation
class Simulation
{
public:
    Simulation() = delete;
    Simulation(const char& boundary_type, const SimulationCellGrid& initial_grid, const int& temperature);
    void step();
    void print_grid(std::ostream& stream) const;

private:
    SimulationCellGrid grid;
    char boundary_type;
    int time;
    int temperature;
    std::vector<Event> event_list;
    void populate_event_list();
    double calculate_rate(const Event&);
    double calculate_repulsion_energy_change(const Event&);
};

#endif
