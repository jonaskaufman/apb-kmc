import numpy as np


class PixelGrid:
    """ A grid of scalar pixel values, periodic in x and y """

    def __init__(self, width, height):
        self.grid = np.zeros((height, width))
        self.width = width
        self.height = height

    def get_value(self, x, y):
        return self.grid[y % self.height][x % self.width]

    def set_value(self, x, y, value):
        self.grid[y % self.height][x % self.width] = value

    def set_grid(self, grid):
        grid = np.array(grid)
        height, width = grid.shape
        assert self.width == width
        assert self.height == height
        self.grid = grid

    def get_smooth_grid(self, sigma):
        window = 4*sigma
        smooth_grid = PixelGrid(self.width, self.height)
        for y in range(self.height):
            for x in range(self.width):
                values = [self.get_value(x, y+dy)*gaussian(0, sigma, dy)
                          for dy in range(-window, window+1)]
                smooth_grid.set_value(x, y, sum(values))
        return smooth_grid

    def get_horizontal_averages(self):
        horizontal_averages = []
        for y in range(self.height):
            horizontal_averages.append(np.mean(self.grid[y]))
        return horizontal_averages


def gaussian(mu, sigma, x):
    return np.exp(-((x-mu)**2)/(2*(sigma**2)))/(sigma*np.sqrt(2*np.pi))

def periodic_smooth(profile, sigma):
    smooth_profile = []
    window = 4*sigma
    length = len(profile)
    for i in range(length):
        values = [profile[(i+j)%length]*gaussian(0, sigma, j) for j in range(-window, window+1)]
        smooth_profile.append(sum(values))
    return smooth_profile

def get_spacing_pixel_grid(phase_grid):
    spacing_grid = PixelGrid(phase_grid.width, phase_grid.height)
    for y in range(phase_grid.height):
        for x in range(phase_grid.width):
            phase = phase_grid.get_value(x, y)
            for dy_up in range(phase_grid.height):
                if phase != phase_grid.get_value(x, y+dy_up):
                    break
            for dy_down in range(phase_grid.height):
                if phase != phase_grid.get_value(x, y-dy_down):
                    break
            k = dy_up+dy_down-1
            spacing_grid.set_value(x, y, k)
    return spacing_grid


def get_composition_pixel_grid(spacing_grid, boundary_type):
    composition_grid = PixelGrid(spacing_grid.width, spacing_grid.height)
    for y in range(spacing_grid.height):
        for x in range(spacing_grid.width):
            k = spacing_grid.get_value(x, y)
            if boundary_type == '-':
                composition = k/(2*k+1)
            elif boundary_type == '+':
                composition = k/(2*k-1)
            composition_grid.set_value(x, y, composition)
    return composition_grid
