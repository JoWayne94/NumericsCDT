"""
Inviscid Burgers' equation types

\partial u / \partial t + div(F) = 0

where the flux, F, is

F = 0.5 * u^2,

and u is a scalar.

Author: JWT
"""

def burgers_flux(a, u):

    return 0.5 * u * u


def inviscid_burgers_1d(dt, velocity_vector, numerical_flux, grid):
    """
    @:brief Inviscid Burgers' flux computations

    :param dt:              Time-step size
    :param grid:            Grid variable object
    :param numerical_flux:  Numerical flux of choice
    :param velocity_vector: Constant in time driving velocity field
    :return: 0.5 * u * u
    """

    # numerical_flux(grid, velocity_vector * 0.5 * grid.variables_data.values)

    grid.variables_data.x_fluxes = numerical_flux(grid.variables_data.values,
                                                  burgers_flux(velocity_vector, grid.variables_data.values),
                                                  grid.dx, dt)


def inviscid_burgers_2d(dt, velocity_vector, numerical_flux, grid):
    """
    @:brief Scalar inviscid Burgers' flux computations

    :param dt:              Time-step size
    :param grid:            Grid variable object
    :param numerical_flux:  Numerical flux of choice
    :param velocity_vector: Velocity field that drives the flow
    :return: 0.5 * u * u
    """

    # Loop through rows
    for j in range(grid.ny):

        grid.variables_data.x_fluxes[j, :] = numerical_flux(grid.variables_data.values[j + 1, :],
                                                            burgers_flux(velocity_vector[j + 1, :, 0], grid.variables_data.values[j + 1, :]),
                                                            grid.dx, dt)

    # Loop through columns
    for i in range(grid.nx):

        grid.variables_data.y_fluxes[:, i] = numerical_flux(grid.variables_data.values[:, i + 1],
                                                            burgers_flux(velocity_vector[:, i + 1, 1], grid.variables_data.values[:, i + 1]),
                                                            grid.dy, dt)
