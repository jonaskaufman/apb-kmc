#include "wrapper.hpp"
#include "grid.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <fstream>
#include <numeric>
#include <vector>

std::vector<int> spacings_for_sinusoidal_composition(BOUNDARY_TYPE boundary_type,
                                                     double composition_average,
                                                     double composition_amplitude,
                                                     double target_physical_height)
{
    // Convert the desired average and amplitude of composition to boundary spacing values
    double spacing_average = average_spacing_from_composition(boundary_type, composition_average);
    double spacing_amplitude =
        std::abs((average_spacing_from_composition(boundary_type, composition_average + composition_amplitude) -
                  average_spacing_from_composition(boundary_type, composition_average - composition_amplitude))) /
        2;
    // Make sure boundary spacing will never be too small
    if (spacing_average - spacing_amplitude < 1)
    {
        std::cerr << "Amplitude too large!" << std::endl;
        return std::vector<int>();
    }

    // Accounting for the added physical height of boundaries themselves,
    // calculate the approximate number of boundaries and total cell height
    double additional_boundary_height = (boundary_type == BOUNDARY_TYPE::MINUS) ? 0.5 : -0.25;
    double spacing_average_adjusted = spacing_average + additional_boundary_height;
    int n_boundaries = std::round((double)target_physical_height / spacing_average_adjusted);
    if (n_boundaries % 2 != 0) // Number of boundaries must be even to ensure periodicity
    {
        n_boundaries += 1;
    }
    int cell_height = std::round(n_boundaries * spacing_average);
    double physical_height = cell_height + additional_boundary_height * n_boundaries;

    // Obtain raw boundary spacing values by sampling sinusoidal composition profile
    std::vector<double> raw_spacings(n_boundaries, -1);
    double phase_shift = (boundary_type == BOUNDARY_TYPE::MINUS) ? 0.0 : M_PI;
    double y_position = 0.0;
    for (int i = 0; i < n_boundaries; i++)
    {
        double t = 2 * M_PI * y_position / physical_height;

        double composition_sample = composition_amplitude * std::cos(t + phase_shift) + composition_average;
        raw_spacings[i] = average_spacing_from_composition(boundary_type, composition_sample);
        if (i < n_boundaries / 2)
        {
            y_position += (raw_spacings[i] + additional_boundary_height);
        }
        else
        {
            y_position -= (raw_spacings[i] + additional_boundary_height);
        }
    }

    // Scale the raw spacing values to get closer to the desired cell height
    double raw_spacings_total = std::accumulate(raw_spacings.begin(), raw_spacings.end(), 0.0);
    double spacing_scaling = (double)cell_height / raw_spacings_total;
    std::vector<int> spacings(n_boundaries, -1);
    for (int i = 0; i < n_boundaries; i++)
    {
        spacings[i] = std::round(spacing_scaling * raw_spacings[i]); // Final spacing values must be integers
    }

    // Make sure cell height is even to ensure periodicity
    cell_height = std::accumulate(spacings.begin(), spacings.end(), 0);
    if (cell_height % 2 != 0)
    {
        spacings[0] += 1;
    }
    return spacings;
}

double average_spacing_from_composition(BOUNDARY_TYPE boundary_type, double composition)
{
    if (boundary_type == BOUNDARY_TYPE::MINUS)
    {
        return (composition / (1 - 2 * composition));
    }
    else
    {
        return (composition / (4 * composition - 2));
    }
}

double physical_height(BOUNDARY_TYPE boundary_type, const std::vector<int>& boundary_cell_spacings)
{
    double cell_height = std::accumulate(boundary_cell_spacings.begin(), boundary_cell_spacings.end(), 0);
    return cell_height + ((boundary_type == BOUNDARY_TYPE::MINUS) ? 0.5 : -0.25) * boundary_cell_spacings.size();
}

void print_initialization_report(BOUNDARY_TYPE boundary_type,
                                 const std::vector<int>& initial_spacings,
                                 double target_physical_height,
                                 double target_composition,
                                 std::ostream& output_stream)
{
    output_stream << "Boundary type: " << (boundary_type == BOUNDARY_TYPE::MINUS ? "-" : "+") << std::endl;
    output_stream << "Initial spacings: ";
    for (auto& a : initial_spacings)
    {
        output_stream << a << " ";
    }
    output_stream << std::endl;
    int cell_height = std::accumulate(initial_spacings.begin(), initial_spacings.end(), 0);
    double true_height = physical_height(boundary_type, initial_spacings);
    output_stream << "Cell height: " << cell_height << std::endl;
    output_stream << "True height (including boundaries): " << true_height << " (Target: " << target_physical_height
                  << ")" << std::endl;
    int n_boundaries = initial_spacings.size();
    double spacing_average = (double)cell_height / (double)n_boundaries;
    double target_spacing_average = average_spacing_from_composition(boundary_type, target_composition);
    output_stream << "Average cell spacing: " << spacing_average << " (Target: " << target_spacing_average << ")"
                  << std::endl;
    return;
}

SimulationWrapper::SimulationWrapper(BOUNDARY_TYPE boundary_type,
                                     int width,
                                     const std::vector<int>& initial_spacings,
                                     double temperature)
    : boundary_type{boundary_type}, width{width}, initial_spacings{initial_spacings}, temperature{temperature}
{
}

void SimulationWrapper::perform_single(int total_passes,
                                       int print_interval,
                                       std::ofstream& phase_grid_file_stream,
                                       std::ofstream& composition_grid_file_stream,
                                       std::ofstream& composition_profile_file_stream)
{
    Simulation simulation = setup();
    for (int n = 0; n < total_passes; n++)
    {
        if (n % print_interval == 0)
        {
            phase_grid_file_stream << simulation.get_time() << std::endl;
            composition_grid_file_stream << simulation.get_time() << std::endl;
            composition_profile_file_stream << simulation.get_time() << std::endl;
            print_pixel_grid(simulation.get_phase_pixel_grid(), phase_grid_file_stream);
            print_pixel_grid(simulation.get_composition_pixel_grid(), composition_grid_file_stream);
            for (auto& a : simulation.get_average_composition_profile())
            {
                composition_profile_file_stream << a << " ";
            }
            composition_profile_file_stream << std::endl;
        }
        simulation.pass();
    }
}

Simulation SimulationWrapper::setup() const
{
    SimulationCellGrid grid(width, initial_spacings, stagger_grid());
    Simulation simulation(boundary_type, grid, temperature);
    return simulation;
}

bool SimulationWrapper::stagger_grid() const { return (boundary_type == BOUNDARY_TYPE::MINUS); }
