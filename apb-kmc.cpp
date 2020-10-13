#include "wrapper.hpp"
#include <fstream>
#include <string>

int main(int argc, char* argv[])
{
    if (argc < 10)
    {
        std::cerr << "Not enough arguments provided. Arguments are (in order):" << std::endl;
        std::cerr << "  boundary type ('-'  or '+')" << std::endl;
        std::cerr << "  grid width (integer)" << std::endl;
        std::cerr << "  grid height (integer)" << std::endl;
        std::cerr << "  composition average (real)" << std::endl;
        std::cerr << "  composition amplitude (real)" << std::endl;
        std::cerr << "  temperature (real, kelvin)" << std::endl;
        std::cerr << "  number of simulations (integer)" << std::endl;
        std::cerr << "  number of passes (integer)" << std::endl;
        std::cerr << "  print interval (integer)" << std::endl;
        std::cerr << "  " << std::endl;
        return 1;
    }
    BOUNDARY_TYPE boundary_type = (*argv[1] == '-' ? BOUNDARY_TYPE::MINUS : BOUNDARY_TYPE::PLUS);
    int width = std::atoi(argv[2]);
    int height = std::atoi(argv[3]);
    double composition_average = std::stod(argv[4]);
    double composition_amplitude = std::stod(argv[5]);
    double temperature = std::stod(argv[6]);
    int total_simulations = std::atoi(argv[7]);
    int total_passes = std::atoi(argv[8]);
    int print_interval = std::atoi(argv[9]);
    std::cout << "done" << std::endl;

    std::vector<int> initial_spacings =
        spacings_for_sinusoidal_composition(boundary_type, composition_average, composition_amplitude, height);

    SimulationWrapper wrapper(boundary_type, width, initial_spacings, temperature);
    std::ofstream composition_profile_file;
    composition_profile_file.open("composition_profile.out");
    if (total_simulations == 1)
    {
        std::ofstream phase_grid_file;
        phase_grid_file.open("phase_grid.out");
        std::ofstream composition_grid_file;
        composition_profile_file.open("composition_grid.out");
        wrapper.perform_single(total_passes, print_interval, phase_grid_file, composition_grid_file,
                               composition_profile_file);
        phase_grid_file.close();
        composition_grid_file.close();
    }
    else
    {
        wrapper.perform_set(total_simulations, total_passes, print_interval, composition_profile_file);
    }
    composition_profile_file.close();
    return 0;
};
