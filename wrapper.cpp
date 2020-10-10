#include "wrapper.hpp"
#include <cmath>
#include <vector>

std::vector<int> spacings_for_sinusoidal_composition(const BOUNDARY_TYPE& boundary_type,
                                                     const double& composition_average,
                                                     const double& composition_amplitude,
                                                     const int& target_height)
{
    std::vector<int> spacing_samples(target_height, 0);
    for (int y = 0; y < target_height; y++)
    {
        double composition = composition_amplitude * std::cos(y * 2 * M_PI / target_height) + composition_average;
        int spacing = 0;
        if (boundary_type == BOUNDARY_TYPE::MINUS)
        {
            spacing = std::round(composition / (1 - 2 * composition));
        }
        else if (boundary_type == BOUNDARY_TYPE::PLUS)
        {
            spacing = std::round(composition / (2 * composition - 1));
        }
        spacing_samples[y] = spacing;
    }
    std::vector<int> spacings;
    int y = 0;
    while (y < target_height / 2)
    {
        int next_spacing = spacing_samples[y];
        spacings.push_back(next_spacing);
        y += next_spacing;
    }
    int length = spacings.size();
    for (int i = 0; i < length; i++)
    {
        spacings.push_back(spacings[length - 1 - i]);
    }
    return spacings;
}

SimulationWrapper::SimulationWrapper(const BOUNDARY_TYPE& boundary_type,
                                     const int& width,
                                     const std::vector<int>& initial_spacings,
                                     const double& temperature)
    : boundary_type{boundary_type}, width{width}, initial_spacings{initial_spacings}, temperature{temperature}
{
}

void SimulationWrapper::perform_single(const int& total_steps,
                                       const int& print_interval,
                                       std::ofstream& output_file_stream)
{
    Simulation simulation = setup();
    for (int n = 0; n < total_steps; n++)
    {
        if (n % print_interval == 0)
        {
            output_file_stream << simulation.get_time() << std::endl;

//            simulation.print_phase_pixel_grid(output_file_stream);
        }
        simulation.step();
    }
}

void SimulationWrapper::perform_set(const int& total_simulations,
                                    const int& total_steps,
                                    const int& print_interval,
                                    std::ofstream& output_file_stream)
{
    std::vector<double> times;
    std::vector<std::vector<std::vector<double>>> average_compositions;
    for (int k = 0; k < total_simulations; k++)
    {
        std::cout << "Running simulation " << k << "..." << std::endl;
        // TODO could initialize this vector based on known size
        average_compositions.emplace_back(std::vector<std::vector<double>>());
        Simulation simulation = setup();
        for (int n = 0; n < total_steps; n++)
        {
            if (n % print_interval == 0)
            {
                if (k == 0)
                {
                    times.emplace_back(simulation.get_time());
                }
                average_compositions[k].emplace_back(simulation.average_horizontal_composition_pixels());
            }
            simulation.step();
        }
        std::cout << "Done" << std::endl;
    }
    std::cout << "Writing results file..." << std::endl;
    for (int t = 0; t < times.size(); t++)
    {
        output_file_stream << times[t] << std::endl;
        for (int k = 0; k < total_simulations; k++)
        {
            for (auto& a : average_compositions[k][t])
            {
                output_file_stream << a << " ";
            }
            output_file_stream << "\n";
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
