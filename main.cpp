#include "simulation.hpp"
#include <fstream>
#include <string>

int main()
{
    int width = 5;
    std::vector<int> spacings{12, 11, 9, 7, 6, 6, 5, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3,  3,
                              3,  3,  3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 6, 6, 7, 9, 11, 12};
    // std::vector<int> spacings{6, 5, 4, 4, 3, 3, 4, 4, 5, 6};

    SimulationCellGrid grid(width, spacings, true);
    // grid.print_pixel_grid(std::cout);

    int steps = 10000000;
    int print_interval = 100000;
    int temperature = 300;

    int n_simulations = 100;

    std::vector<double> times;
    std::vector<std::vector<std::vector<double>>> average_spacings;
    for (int k = 0; k < n_simulations; k++)
    {
        average_spacings.emplace_back(std::vector<std::vector<double>>());
        Simulation simulation(BOUNDARY_TYPE::MINUS, grid, temperature);
        for (int n = 0; n < steps; n++)
        {
            if (n % print_interval == 0)
            {
                if (k == 0)
                {
                    times.emplace_back(simulation.get_time());
                }
                average_spacings[k].emplace_back(simulation.get_horizontal_pixel_average_spacings());
            }
            simulation.step();
        }
    }

    for (int t = 0; t < times.size(); t++)
    {
        std::cout << times[t] << std::endl;
        for (int k = 0; k < n_simulations; k++)
        {
            for (auto& a : average_spacings[k][t])
            {
                std::cout << a << " ";
            }
            std::cout << std::endl;
        }
    }
}
