import numpy as np


def parse_profile_file(file_name):
    """ Parse a file containing a set of profiles at each time stamp """
    times = []
    profile_sets = []
    with open(file_name, 'r') as f:
        profile_set = []
        while True:
            split_line = f.readline().split()
            if len(split_line) > 1:
                profile_set.append(list(map(float, split_line)))
            else:
                if profile_set:
                    profile_sets.append(profile_set)
                if len(split_line) == 1:
                    times.append(float(split_line[0]))
                    profile_set = []
                else:   # end of file
                    break
    return times, profile_sets


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
