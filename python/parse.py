import sys
import numpy as np
from grid import *


def parse_results(file_name):
    times = []
    phase_grids = []
    with open(file_name, 'r') as f:
        next_grid = []
        while True:
            split_line = f.readline().split()
            if len(split_line) > 1:
                values = list(map(int, split_line))
                next_grid.append(values)
            else:
                if next_grid:
                    phase_grids.append(
                        PixelGrid(len(next_grid[0]), len(next_grid)))
                    phase_grids[-1].set_grid(next_grid)
                if len(split_line) == 1:
                    times.append(float(split_line[0]))
                    next_grid = []
                else:   # end of file
                    break
    return times, phase_grids


file_name = sys.argv[1]
times, phase_grids = parse_results(file_name)
