#include "wrapper.hpp"
#include <fstream>
#include <string>

int main()
{
    BOUNDARY_TYPE boundary_type = BOUNDARY_TYPE::MINUS;
    double composition_average = 5.0/11.0;
    double composition_amplitude = 0.025;
        int target_height = 200;
    std::vector<int> initial_spacings = spacings_for_sinusoidal_composition(boundary_type, composition_average, composition_amplitude, target_height);
int width = 20;
double temperature = 300;
SimulationWrapper wrapper(boundary_type, width, initial_spacings, temperature);

std::string output_file_name = "results.out";
std::ofstream output_file;
output_file.open(output_file_name);
wrapper.perform_set(100, 10000000, 100000, output_file);
output_file.close();
}
