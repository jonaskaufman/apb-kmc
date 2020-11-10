import numpy as np

grid_x_scaling = 2  # width of one simulation cell, in pixels
grid_y_scaling = 4  # height of one simulation cell, in pixels


def gaussian(mu, sigma, x):
    """ Gaussian function evaluated at x, for average mu and standard deviation sigma """
    return np.exp(-((x-mu)**2)/(2*(sigma**2)))/(sigma*np.sqrt(2*np.pi))


def periodic_gaussian_matrix(dimension, sigma):
    """ Matrix to apply periodic gaussian smoothing to a 1-D profile of given dimension """
    matrix = np.zeros((dimension, dimension))
    window = 4*sigma    # should be enough sigmas
    for i in range(dimension):
        for w in range(-window, window+1):
            j = (i+w) % dimension
            matrix[i, j] += gaussian(0, sigma, w)
    return matrix

