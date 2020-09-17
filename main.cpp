#include "simulation.hpp"


int main()
{
    int width = 5;
    std::vector<int> spacings{12, 11, 9, 7, 6, 6, 5, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3,  3,
                              3,  3,  3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 6, 6, 7, 9, 11, 12};

SimulationCellGrid grid(width, spacings, true);
grid.print_pixel_grid(std::cout);

int steps = 10000;
int temperature = 300;

}

