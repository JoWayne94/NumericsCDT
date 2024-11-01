"""
Scalar transport equation types

\partial u / \partial t + div(F) = 0

where the flux, F, is

F = u \dot v

and v is the driving velocity vector.

Author: JWT
"""

def transport_flux(a, u):

    return a * u


def scalar_transport_1d(dt, velocity_vector, numerical_flux, grid):
    """
    @:brief Scalar transport flux computations

    :param dt:              Time-step size
    :param grid:            Grid variable object
    :param numerical_flux:  Numerical flux of choice
    :param velocity_vector: Velocity field that drives the flow
    :return: v * u
    """

    # fluxes = numerical_flux(fluxes, values, n, velocity_vector)
    #
    # return fluxes

    # First input should be grid.variables_data.values, but due to upwind conditions it is temporarily
    # written as velocity_vector, should be fixed in the future
    grid.variables_data.x_fluxes = numerical_flux(velocity_vector,
                                                  transport_flux(velocity_vector, grid.variables_data.values),
                                                  grid.dx, dt)


def scalar_transport_2d(dt, velocity_vector, numerical_flux, grid):
    """
    @:brief Scalar transport flux computations

    :param dt:              Time-step size
    :param grid:            Grid variable object
    :param numerical_flux:  Numerical flux of choice
    :param velocity_vector: Velocity field that drives the flow
    :return: v * u
    """

    # Loop through rows
    for j in range(grid.ny):

        grid.variables_data.x_fluxes[j, :] = numerical_flux(velocity_vector[j + 1, :, 0],
                                                            transport_flux(velocity_vector[j + 1, :, 0], grid.variables_data.values[j + 1, :]),
                                                            grid.dx, dt)

    # Loop through columns
    for i in range(grid.nx):

        grid.variables_data.y_fluxes[:, i] = numerical_flux(velocity_vector[:, i + 1, 1],
                                                            transport_flux(velocity_vector[:, i + 1, 1], grid.variables_data.values[:, i + 1]),
                                                            grid.dy, dt)
