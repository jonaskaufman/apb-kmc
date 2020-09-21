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
