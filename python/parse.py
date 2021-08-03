import numpy as np


def parse_profile_file(file_name):
    """ Parse a file containing a single 1-D profile at each time stamp """
    times = []
    profiles = []
    with open(file_name, 'r') as f:
        while True:
            time_line = f.readline()
            if time_line:
                times.append(float(time_line.split()[0]))
                profile_line = f.readline()
                profiles.append(list(map(float, profile_line.split())))
            else:   # end of file
                break
    return times, profiles


def parse_grid_file(file_name):
    """ Parse a file containing a single grid at each time stamp """
    times = []
    grids = []
    with open(file_name, 'r') as f:
        next_grid = []
        while True:
            split_line = f.readline().split()
            if len(split_line) > 1:
                values = list(map(float, split_line))
                next_grid.append(values)
            else:
                if next_grid:
                    grids.append(np.array(next_grid))
                if len(split_line) == 1:
                    times.append(float(split_line[0]))
                    next_grid = []
                else:   # end of file
                    break
    return times, grids
