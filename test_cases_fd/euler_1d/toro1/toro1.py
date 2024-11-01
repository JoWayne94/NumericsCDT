"""
Compressible Euler equation set, Toro test 1

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
from src.solvers.compressibleeuler import *


if __name__ == '__main__':
    """
    main()
    """

    if NO_OF_DIMENSIONS == 1:
        grid = Grid1D

        if SPATIAL_SCHEME == 'Finite difference':
            spatial_scheme = central_difference_1d
        # elif SPATIAL_SCHEME == 'CABARET':
        # spatial_scheme = CABARET
        else:
            raise NotImplementedError('Spatial discretisation not implemented.')

        if EQUATION_TYPE == 'Scalar transport':
            eqn_type = scalar_transport_1d
        elif EQUATION_TYPE == 'Inviscid Burgers':
            eqn_type = inviscid_burgers_1d
        elif EQUATION_TYPE == 'Compressible Euler':
            eqn_type = compressible_euler_1d
        else:
            raise NotImplementedError('Equation type is not implemented.')

    elif NO_OF_DIMENSIONS == 2:
        grid = Grid2D

        if SPATIAL_SCHEME == 'Finite difference':
            spatial_scheme = central_difference_2d
        # elif SPATIAL_SCHEME == 'CABARET':
        # spatial_scheme = CABARET
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


    """ Initialise conservative variables """
    u_1 = grid(NO_OF_GRID_POINTS, DOMAIN_BOUNDARIES, BOUNDARY_CONDITIONS_TYPE)
    u_2 = grid(NO_OF_GRID_POINTS, DOMAIN_BOUNDARIES, BOUNDARY_CONDITIONS_TYPE)
    u_3 = grid(NO_OF_GRID_POINTS, DOMAIN_BOUNDARIES, BOUNDARY_CONDITIONS_TYPE)

    '''======PLEASE INPUT YOUR INITIAL CONDITIONS BELOW======'''
    # x coordinates
    x = u_1.points_data.coordinates
    # Initial conditions defined with a numpy array
    rho_ic = np.where(x <= 0.5, 1., 0.125)
    u_ic = np.where(x <= 0.5, 0., 0.)
    p_ic = np.where(x <= 0.5, 1., 0.1)
    '''======================================================'''

    u_1.variables_data.values, u_2.variables_data.values, u_3.variables_data.values = prim_to_cons(rho_ic, u_ic, p_ic)
    v = np.ones(u_1.nx + 2)

    """ Compute speed of sound and calculate time-step constraint """
    c = speed_of_sound(p_ic, rho_ic)
    dt = CFL * u_1.dx / np.max(np.abs(u_ic) + c)

    """ Plotting parameters and visualisations """
    lw = 1.5
    y_bottom = 0.
    y_top = 1.

    # Initial conditions
    plt.plot(x[1:-1], u_1.variables_data.values[1:-1], 'k', label='IC')
    plt.legend(loc='best')
    plt.ylabel(r'$\rho$')
    # plt.axhline(H, linestyle=':', color='black')
    plt.ylim([y_bottom, y_top])
    plt.pause(1)

    t = 0 # Initial time
    while abs(FINAL_TIME - t) > 1.e-7:

        if t + dt >= FINAL_TIME: dt = FINAL_TIME - t

        """ Equation type """
        eqn_type(dt, v, numerical_flux, u_1, u_2, u_3)

        """ Spatial discretisation """
        spatial_scheme(u_1, 1.)
        spatial_scheme(u_2, 1.)
        spatial_scheme(u_3, 1.)

        """ Time evolution """
        forward_euler(u_1, dt)
        forward_euler(u_2, dt)
        forward_euler(u_3, dt)


        # Replot
        plt.cla()
        plt.plot(x[1:-1], u_1.variables_data.values[1:-1], 'b', label='Time = ' + '%.2f' % (t + dt) + ' s', linewidth=lw)
        plt.title(r'$\sigma = $' + str(CFL))
        plt.legend(loc='lower left')
        plt.xlabel('x')
        plt.ylabel(r'$\rho$')
        plt.ylim([y_bottom, y_top])
        plt.pause(0.01)

        t += dt

        rho_ic, u_ic, p_ic = cons_to_prim(u_1.variables_data.values, u_2.variables_data.values, u_3.variables_data.values)
        c = speed_of_sound(p_ic, rho_ic)
        dt = CFL * u_1.dx / np.max(np.abs(u_ic) + c)

    """ Analytical solution """
    # ac = np.zeros_like(x)
    # for i in range(x.shape[0]):
    #
    #     if 0. <= x[i] < FINAL_TIME:
    #         ac[i] = float(x[i]) / FINAL_TIME
    #     if FINAL_TIME <= x[i] < 0.5 * FINAL_TIME + 1.:
    #         ac[i] = 1.
    #
    # plt.plot(x[1:-1], ac[1:-1], 'k--', label='Exact', linewidth=lw)
    # plt.legend(loc='best')
    # plt.ylabel(r'$u$')
    # plt.title(r'$\sigma = $' + str(CFL) + ', ' + r'$T = L^{-1} u$')
    # plt.axhline(0, linestyle=':', color='black')
    # plt.grid()
    plt.show()
