#include "simulation.hpp"
#include <random>
#include <stdexcept>

// TODO check order is same as declaration
Simulation::Simulation(BOUNDARY_TYPE boundary_type, const SimulationCellGrid& initial_grid, double temperature)
    : boundary_minus{boundary_type == BOUNDARY_TYPE::MINUS},
      grid_ptr{std::make_shared<SimulationCellGrid>(initial_grid)},
      temperature{temperature},
      time{0}
{
    // Check if grid is compatible with boundary type
    bool correct_staggered = boundary_minus;
    if (grid_ptr->staggered != correct_staggered)
    {
        throw std::runtime_error("Boundary type not compatible with grid.");
    }

    // Initialize event list
    event_list_ptr = std::make_shared<const std::vector<Event>>(generate_event_list());

    // Initialize rate calculator
    rate_calculator_ptr = std::make_shared<EventRateCalculator>(boundary_type, temperature, event_list_ptr, grid_ptr);

    // Initialize event selector
    double rate_upper_bound = 1.0;
    selector_ptr = std::make_unique<lotto::RejectionEventSelector<ID, EventRateCalculator>>(
        rate_calculator_ptr, rate_upper_bound, generate_event_id_list());
}

void Simulation::step()
{
    // Select event
    auto event_id_and_time = selector_ptr->select_event();
    ID accepted_event_id = event_id_and_time.first;
    double time_step = event_id_and_time.second;

    // Update grid
    for (auto& indices : event_list_ptr->at(accepted_event_id))
    {
        grid_ptr->flip_cell_phase(indices.first, indices.second);
    }

    // Update time step
    if (time_step < 0)
    {
        throw std::runtime_error("Time step is invalid.");
    }
    time += time_step;
}

void Simulation::pass()
{
    for (int i = 0; i < event_list_ptr->size(); ++i)
    {
        step();
    }
}

PixelGrid Simulation::get_phase_pixel_grid() const
{
    const PixelGrid phase_grid = grid_ptr->get_phase_pixel_grid();
    int width = phase_grid.size();
    int height = phase_grid[0].size();
    int scaled_height = 4 * height;
    PixelGrid scaled_phase_grid(width, std::vector<double>(scaled_height, 0));
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            for (int dy = 0; dy < 4; dy++)
            {
                scaled_phase_grid[x][4 * y + dy] = phase_grid[x][y];
            }
        }
    }
    return scaled_phase_grid;
}

PixelGrid Simulation::get_composition_pixel_grid() const
{
    const PixelGrid phase_grid = grid_ptr->get_phase_pixel_grid();
    int width = phase_grid.size();
    int phase_grid_height = phase_grid[0].size();
    // Find the y coordinate of each boundary at x = 0
    std::vector<int> boundary_y_origins;
    for (int y = 0; y < phase_grid_height; y++)
    {
        double current_phase = phase_grid[0][y];
        double previous_phase = phase_grid[0][modulo(y - 1, phase_grid_height)];
        if (current_phase != previous_phase)
        {
            boundary_y_origins.push_back(y);
        }
    }

    // zeta plus currently wrong!
    int n_boundaries = boundary_y_origins.size();
    int composition_grid_height =
        (boundary_minus ? (4 * phase_grid_height + 2 * n_boundaries) : (4 * phase_grid_height - n_boundaries));
    PixelGrid composition_grid(width, std::vector<double>(composition_grid_height, 0.5));
    // Follow each boundary along the x direction and update composition accordingly
    for (int b = 0; b < n_boundaries; b++)
    {
        int y_phase = boundary_y_origins[b];
        int y_comp = (boundary_minus ? (4 * y_phase + 2 * b) : (4 * y_phase - b));
        for (int x = 0; x < width; x++)
        {
            if (boundary_minus)
            {
                composition_grid[x][modulo(y_comp, composition_grid_height)] = 0.0;
                composition_grid[x][modulo(y_comp + 1, composition_grid_height)] = 0.0;
            }
            else
            {
                composition_grid[x][modulo(y_comp, composition_grid_height)] = 1.0;
            }
            int x_next = modulo(x + 1, width);
            int y_phase_previous = modulo(y_phase - 1, phase_grid_height);
            int y_phase_current = modulo(y_phase, phase_grid_height);
            if (phase_grid[x][y_phase_current] != phase_grid[x_next][y_phase_current])
            {
                y_phase += 1;
                y_comp += 4;
            }
            else if (phase_grid[x][y_phase_previous] != phase_grid[x_next][y_phase_previous])
            {
                y_phase -= 1;
                y_comp -= 4;
            }
        }
    }
    return composition_grid;
}

std::vector<double> Simulation::get_average_composition_profile() const
{
    return average_pixels_horizontally(get_composition_pixel_grid());
}

std::vector<Event> Simulation::generate_event_list() const
{
    std::vector<Event> event_list;

    for (int y = 0; y < grid_ptr->height; ++y)
    {
        for (int x = 0; x < grid_ptr->width; ++x)
        {
            event_list.push_back({std::make_pair(x, y)});
            if (!boundary_minus)
            {
                event_list.push_back({std::make_pair(x, y), std::make_pair(modulo(x + 1, grid_ptr->width), y)});
            }
        }
    }
    return event_list;
}

std::vector<ID> Simulation::generate_event_id_list() const
{
    std::vector<ID> event_id_list;
    for (int i = 0; i < event_list_ptr->size(); ++i)
    {
        event_id_list.push_back(i);
    }
    return event_id_list;
}

std::map<Coordinates, std::vector<ID>> Simulation::generate_coordinates_ids_map() const
{
    std::map<Coordinates, std::vector<ID>> coordinates_ids_map;
    for (const ID& event_id : generate_event_id_list())
    {
        for (const Coordinates& coordinates : event_list_ptr->at(event_id))
        {
            coordinates_ids_map[coordinates].push_back(event_id);
        }
    }
    return coordinates_ids_map;
}

std::set<Coordinates> Simulation::generate_impact_neighborhood(const Event& impacted_event) const
{
    std::set<Coordinates> neighborhood;
    for (const Coordinates& impacted_coordinates : impacted_event)
    {
        neighborhood.insert(impacted_coordinates);

        // zeta minus
        //          NN
        //      NW      NE
        //  W       0       E
        //      SW      SE
        //          SS

        // zeta plus
        //          NN
        //      NW  N   NE
        //  WW  W   0   E   EE
        //      SW  S   SE
        //          SS
    }
    return neighborhood;
}

std::map<ID, std::vector<ID>> Simulation::generate_impact_table() const
{
    std::vector<ID> event_id_list = generate_event_id_list();
    std::map<Coordinates, std::vector<ID>> coordinates_ids_map = generate_coordinates_ids_map();
    std::map<ID, std::vector<ID>> impact_table;
    // Go through impacted events
    for (const ID& impacted_event_id : event_id_list)
    {
        // Go through cells that impact event's rate
        for (const Coordinates& impacter_coordinates :
             generate_impact_neighborhood(event_list_ptr->at(impacted_event_id)))
        {
            // Go through events that contain those cells
            for (const ID& impacter_event_id : coordinates_ids_map.at(impacter_coordinates))
            {
                // Add impacted event to appropriate places in impact table
                impact_table[impacter_event_id].push_back(impacted_event_id);
            }
        }
    }
    return impact_table;
};

