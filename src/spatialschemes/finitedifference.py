"""
Finite difference spatial schemes - first-order as default

Author: JWT
"""
import numpy as np


# def one_sided_difference(grid, coefficient, direction='forward'):
#
#     first_derivative = -coefficient * np.diff(grid.variables_data.values) / grid.dx
#
#     if direction == 'forward': return first_derivative[1:]
#     elif direction == 'backward': return first_derivative[:-1]
#     else: raise ValueError('Direction must be either forward or backward.')


def conservative_fd_1d(grid, coefficient):

    grid.variables_data.x_first_derivative = -coefficient * (grid.variables_data.x_fluxes[1:] - grid.variables_data.x_fluxes[:-1]) / grid.dx


def conservative_fd_2d(grid, coefficient):

    # first_derivative = -coefficient * (fluxes[1:] - fluxes[:-1]) / delta
    grid.variables_data.x_first_derivative = -coefficient * (grid.variables_data.x_fluxes[:, 1:] - grid.variables_data.x_fluxes[:, :-1]) / grid.dx
    grid.variables_data.y_first_derivative = -coefficient * (grid.variables_data.y_fluxes[1:, :] - grid.variables_data.y_fluxes[:-1, :]) / grid.dy

