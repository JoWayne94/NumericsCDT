"""
A sine wave translating to the left or right

Author: JWT
"""

# Import third-party libraries
import sys, os
import matplotlib.pyplot as plt

# Import user inputs
from setup import *

# Configure system path
if (sys.platform[:3] == 'win') or (sys.platform[:3] == 'Win'):
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..\..')))
else:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

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


if __name__ == '__main__':
    """
    main()
    """

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

    if FLUX_SCHEME == 'Upwind':
        numerical_flux = upwind
    elif FLUX_SCHEME == 'FORCE':
        numerical_flux = force
    else:
        raise NotImplementedError('Numerical flux of choice is not implemented.')


    var = grid(NO_OF_GRID_POINTS, DOMAIN_BOUNDARIES, BOUNDARY_CONDITIONS_TYPE)

    '''======PLEASE INPUT YOUR INITIAL CONDITIONS BELOW======'''
    # x coordinates
    x = var.points_data.coordinates
    # Initial conditions defined with a numpy array
    ic = np.sin(x)
    '''======================================================'''

    var.variables_data.values = ic
    v = -1. * np.ones(var.nx + 2)
    dt = CFL * var.dx / np.max(np.abs(v))

    """ Plotting parameters and visualisations """
    lw = 1.5
    y_bottom = -1.
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
        forward_euler(var, dt)

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
    exact_solution = np.sin(x[1:-1] - v[1:-1] * FINAL_TIME)
    plt.plot(x[1:-1], exact_solution, 'k--', label='Exact soln', linewidth=lw)
    plt.legend(loc='best')
    plt.ylabel(r'$\phi$')
    plt.title(r'$CFL = $' + str(CFL) + ', ' + r'$T = L^{-1} u$')
    plt.axhline(0, linestyle=':', color='black')
    plt.grid()
    plt.show()
