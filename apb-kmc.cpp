#include "submodules/nlohmann-json/include/nlohmann/json.hpp"
#include "wrapper.hpp"
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char* argv[])
{
    // Read input settings from json file
    std::filesystem::path input_file_name = "input.json";
    if (!std::filesystem::exists(input_file_name))
    {
        throw std::runtime_error("Error: input.json file missing");
    }
    nlohmann::json input_settings;
    std::ifstream input_stream(input_file_name);
    input_stream >> input_settings;

    // Parse input settings
    BOUNDARY_TYPE boundary_type = (input_settings["boundary_type"] == "-" ? BOUNDARY_TYPE::MINUS : BOUNDARY_TYPE::PLUS);
    std::vector<double> target_grid_dimensions = input_settings["target_grid_dimensions"];
    int width = target_grid_dimensions[0];
    double target_height = target_grid_dimensions[1];
    double temperature = input_settings["temperature"];
    int total_passes = input_settings["total_passes"];
    int print_interval = input_settings["print_interval"];

    // Get initialization settings and calculate initial boundary spacings
    std::vector<int> initial_spacings;
    auto initialization_settings = input_settings["initialization"];
    double composition_average = initialization_settings["composition_average"];
    std::string initialization_mode = initialization_settings["mode"];
    if (initialization_mode == "sinusoidal")
    {
        double composition_amplitude = initialization_settings["composition_amplitude"];
        initial_spacings = spacings_for_sinusoidal_composition(boundary_type, composition_average,
                                                               composition_amplitude, target_height);
    }
    else
    {
        throw std::runtime_error("Error: initalization mode '" + initialization_mode + "' is not supported");
    }
    print_initialization_report(boundary_type, initial_spacings, target_height, composition_average, std::cout);

    // Set up and run simulation
    try
    {
        SimulationWrapper wrapper(boundary_type, width, initial_spacings, temperature);
        std::ofstream phase_grid_file;
        std::ofstream composition_grid_file;
        std::ofstream composition_profile_file;
        phase_grid_file.open("phase_grid.out");
        composition_grid_file.open("composition_grid.out");
        composition_profile_file.open("composition_profile.out");

        std::time_t start, end;
        std::time(&start);
        std::cout << "Running simulation..." << std::endl;
        wrapper.perform_single(total_passes, print_interval, phase_grid_file, composition_grid_file,
                               composition_profile_file);
        phase_grid_file.close();
        composition_grid_file.close();
        composition_profile_file.close();
        std::cout << "Done. Results files written." << std::endl;
        std::time(&end);
        double elapsed = double(end - start);
        std::cout << "Elapsed time: " << elapsed << " seconds" << std::endl;
    }
    catch (const std::runtime_error& error)
    {
        std::cerr << error.what() << std::endl;
    }
    return 0;
};
