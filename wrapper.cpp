#include "wrapper.hpp"
#include "grid.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <fstream>
#include <numeric>
#include <vector>

std::vector<int> spacings_for_sinusoidal_composition(const BOUNDARY_TYPE& boundary_type,
                                                     const double& composition_average,
                                                     const double& composition_amplitude,
                                                     const int& target_height)
{
    double spacing_average = average_spacing_from_composition(boundary_type, composition_average);
    double spacing_amplitude =
        std::abs((average_spacing_from_composition(boundary_type, composition_average + composition_amplitude) -
                  average_spacing_from_composition(boundary_type, composition_average - composition_amplitude))) /
        2;
    if (spacing_average - spacing_amplitude < 1)
    {
        std::cerr << "Amplitude too large!" << std::endl;
        return std::vector<int>();
    }
    int n_boundaries = std::round((double)target_height / spacing_average);
    if (n_boundaries %2 != 0)
    {
        n_boundaries += 1;  
    }
    int height = std::round(n_boundaries * spacing_average);
    std::vector<double> raw_spacings(n_boundaries, -1);
    for (int i = 0; i < n_boundaries; i++)
    {
        double x = 2 * M_PI * (double)i / (double)n_boundaries;
        double arg_scaling = std::sin(x / 2);
        raw_spacings[i] = spacing_amplitude * std::cos((x + M_PI * arg_scaling) / (1 + arg_scaling)) + spacing_average;
    }
    double raw_spacings_total = std::accumulate(raw_spacings.begin(), raw_spacings.end(), 0.0);
    double spacing_scaling = (double)height / raw_spacings_total;
    std::vector<int> spacings(n_boundaries, -1);
    for (int i = 0; i < n_boundaries; i++)
    {
        spacings[i] = std::round(spacing_scaling * raw_spacings[i]);
    }
    height = std::accumulate(spacings.begin(), spacings.end(), 0);
    if (height %2 != 0)
    {
        spacings[0] += 1; 
    }
    return spacings;
}

double average_spacing_from_composition(const BOUNDARY_TYPE& boundary_type, const double& composition)
{
    if (boundary_type == BOUNDARY_TYPE::MINUS)
    {
        return (composition / (1 - 2 * composition));
    }
    else
    {
        return (composition / (2 * composition - 1));
    }
}

SimulationWrapper::SimulationWrapper(const BOUNDARY_TYPE& boundary_type,
                                     const int& width,
                                     const std::vector<int>& initial_spacings,
                                     const double& temperature)
    : boundary_type{boundary_type}, width{width}, initial_spacings{initial_spacings}, temperature{temperature}
{
}

void SimulationWrapper::perform_single(const int& total_passes,
                                       const int& print_interval,
                                       std::ofstream& phase_grid_file_stream,
                                       std::ofstream& composition_grid_file_stream,
                                       std::ofstream& composition_profile_file_stream)
{
    Simulation simulation = setup();
    std::cout << "Running simulation..." << std::endl;
    for (int n = 0; n < total_passes; n++)
    {
        if (n % print_interval == 0)
        {
            phase_grid_file_stream << simulation.get_time() << std::endl;
            composition_grid_file_stream << simulation.get_time() << std::endl;
            composition_profile_file_stream << simulation.get_time() << std::endl;
            print_pixel_grid(simulation.get_phase_pixel_grid(), phase_grid_file_stream);
            print_pixel_grid(simulation.get_composition_pixel_grid(), composition_grid_file_stream);
            for (auto& a : simulation.average_horizontal_composition_pixels())
            {
                composition_profile_file_stream << a << " ";
            }
            composition_profile_file_stream << std::endl;
        }
        simulation.pass();
    }
        std::cout << "Done. Results files written." << std::endl;
}

void SimulationWrapper::perform_set(const int& total_simulations,
                                    const int& total_passes,
                                    const int& print_interval,
                                    std::ofstream& composition_profile_file_stream)
{
    std::vector<double> times;
    std::vector<std::vector<std::vector<double>>> average_compositions(total_simulations, std::vector<std::vector<double>>());
    for (int k = 0; k < total_simulations; k++)
    {
        std::cout << "Running simulation " << k << "..." << std::endl;
        Simulation simulation = setup();
        for (int n = 0; n < total_passes; n++)
        {
            if (n % print_interval == 0)
            {
                if (k == 0)
                {
                    times.push_back(simulation.get_time());
                }
                average_compositions[k].push_back(simulation.average_horizontal_composition_pixels());
            }
            simulation.pass();
        }
        std::cout << "Done" << std::endl;
    }
    std::cout << "Writing results file..." << std::endl;
    for (int t = 0; t < times.size(); t++)
    {
        composition_profile_file_stream << times[t] << std::endl;
        for (int k = 0; k < total_simulations; k++)
        {
            for (auto& a : average_compositions[k][t])
            {
                composition_profile_file_stream << a << " ";
            }
            composition_profile_file_stream << "\n";
        }
    }
    std::cout << "Done" << std::endl;
}

Simulation SimulationWrapper::setup() const
{
    SimulationCellGrid grid(width, initial_spacings, stagger_grid());
    Simulation simulation(boundary_type, grid, temperature);
    return simulation;
}

bool SimulationWrapper::stagger_grid() const { return (boundary_type == BOUNDARY_TYPE::MINUS); }
