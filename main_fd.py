# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
"""
Main module

Note:

    1. Serves as a template to create other solvers
    2. All grids are non-staggered; first-order conservative finite difference on Arakawa A-Grid/finite volume equivalent
    3. Computational domain only limited to quadrilaterals
    4. Grid spacings are constant for each spatial dimension (\delta x is not a function of x)

Future works:

1. Extension to 2D Euler
2. Compressible Euler and SWE solvers not ready yet

Author: JWT
"""

# Import third-party libraries
import sys, os
import matplotlib.pyplot as plt

# Import user settings
from setup import *

# Import grid classes
from src.grid.grid1D import Grid1D
from src.grid.grid2D import Grid2D

# Import flux definitions & spatial and temporal schemes
from src.spatialschemes.finitedifference import *
from src.temporalschemes.multistage import *
from src.flux.numericalfluxes import *

# Import equation types
from src.solvers.scalartransport import *
from src.solvers.inviscidburgers import *

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press ⌘F8 to toggle the breakpoint.


if __name__ == '__main__':
    """
    main()
    """

    # if len(sys.argv) < 2:
    if TEMPORAL_SCHEME == 'Forward Euler':
        temporal_scheme = forward_euler
    else:
        raise NotImplementedError('Temporal integration method not implemented.')

    if FLUX_SCHEME == 'Upwind':
        numerical_flux = upwind
    elif FLUX_SCHEME == 'FORCE':
        numerical_flux = force
    else:
        raise NotImplementedError('Numerical flux of choice is not implemented.')

    if NO_OF_DIMENSIONS == 1:
        grid = Grid1D

        if SPATIAL_SCHEME == 'Finite difference':
            spatial_scheme = conservative_fd_1d
        else:
            raise NotImplementedError('Spatial discretisation not implemented.')

        if EQUATION_TYPE == 'Scalar transport':
            eqn_type = scalar_transport_1d
        elif EQUATION_TYPE == 'Inviscid Burgers':
            eqn_type = inviscid_burgers_1d
        else:
            raise NotImplementedError('Equation type is not implemented.')

    elif NO_OF_DIMENSIONS == 2:
        grid = Grid2D

        if SPATIAL_SCHEME == 'Finite difference':
            spatial_scheme = conservative_fd_2d
        else:
            raise NotImplementedError('Spatial discretisation not implemented.')

        if EQUATION_TYPE == 'Scalar transport':
            eqn_type = scalar_transport_2d
        elif EQUATION_TYPE == 'Inviscid Burgers':
            eqn_type = inviscid_burgers_2d
        else:
            raise NotImplementedError('Equation type is not implemented.')

    else:
        raise NotImplementedError('Number of spatial dimensions is not supported.')


    var = grid(NO_OF_GRID_POINTS, DOMAIN_BOUNDARIES, BOUNDARY_CONDITIONS_TYPE)

    '''======PLEASE INPUT YOUR INITIAL CONDITIONS BELOW======'''
    center = -0.25
    # x coordinates
    x = var.points_data.coordinates
    # Initial conditions defined with a numpy array
    ic = np.where((x - center) % 2. < 0.5, np.power(np.sin(2. * (x - center) * np.pi), 2), 0.)
    # Prescribed velocity field
    v = 1. * np.ones(var.nx + 2)
    '''======================================================'''

    var.variables_data.values = ic
    dt = CFL * var.dx / np.max(np.abs(v))

    """ Plotting parameters and visualisations """
    lw = 1.5
    y_bottom = 0.
    y_top = 1.

    # Initial conditions
    plt.plot(x[1:-1], var.variables_data.values[1:-1], 'k', label='IC')
    plt.legend(loc='best')
    plt.ylabel(r'$\phi$')
    # plt.axhline(H, linestyle=':', color='black')
    plt.ylim([y_bottom, y_top])
    plt.pause(1)

    t = 0 # Initial time
    while abs(FINAL_TIME - t) > 1.e-7:

        if t + dt >= FINAL_TIME: dt = FINAL_TIME - t

        """ Equation type """
        eqn_type(dt, v, numerical_flux, var)

        """ Spatial discretisation """
        spatial_scheme(var, 1.)

        """ Time evolution """
        temporal_scheme(var, dt)

        # Replot
        plt.cla()
        plt.plot(x[1:-1], var.variables_data.values[1:-1], 'b', label='Time = ' + '%.2f' % (t + dt) + ' s', linewidth=lw)
        # exact = np.where((x - u * (n + 1) * dt) % 1. < 0.5, np.power(np.sin(2. * (x - u * (n + 1) * dt) * np.pi), 2), 0.)
        # plt.plot(x, exact, 'k--', label='Analytical')
        plt.title(r'$CFL = $' + str(CFL))
        plt.legend(loc='lower left')
        plt.xlabel('x')
        plt.ylabel(r'$\phi$')
        plt.ylim([y_bottom, y_top])
        plt.pause(0.01)

        t += dt


    """ Analytical solution """
    exact_solution = np.where((x[1:-1] - center) % 2. < 0.5, np.power(np.sin(2. * (x[1:-1] - center) * np.pi), 2), 0.)
    plt.plot(x[1:-1], exact_solution, 'k--', label='Exact soln', linewidth=lw)
    plt.legend(loc='best')
    plt.ylabel(r'$\phi$')
    plt.title(r'$CFL = $' + str(CFL) + ', ' + r'$T = L^{-1} u$')
    plt.axhline(0, linestyle=':', color='black')
    plt.grid()
    plt.show()
