"""
1D inviscid Burgers plots with a sine wave

Author: JWT
"""

# Import third-party libraries
import sys, os
import matplotlib.pyplot as plt
import numpy as np

# Import user settings
from setup import *

# Configure system path
path = os.path.dirname(__file__)

if (sys.platform[:3] == 'win') or (sys.platform[:3] == 'Win'):
    sys.path.append(os.path.abspath(os.path.join(path, '..\..')))
else:
    sys.path.append(os.path.abspath(os.path.join(path, '../..')))


# Import grid classes
from src.grid.grid1D import Grid1D
from src.grid.grid2D import Grid2D

# Import CABARET methods
from src.CABARET.CABARET import *
from src.CABARET.mesh.mesh1D import *

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
        mesh = Mesh1D

        if SPATIAL_SCHEME == 'Finite difference':
            spatial_scheme = conservative_fd_1d
        else:
            raise NotImplementedError('Spatial discretisation not implemented.')

        if EQUATION_TYPE == 'Scalar transport':
            eqn_type = scalar_transport_1d
            eqn_flux = transport_flux
        elif EQUATION_TYPE == 'Inviscid Burgers':
            eqn_type = inviscid_burgers_1d
            eqn_flux = burgers_flux
        else:
            raise NotImplementedError('Equation type is not implemented.')

    elif NO_OF_DIMENSIONS == 2:
        grid = Grid2D
        mesh = Mesh1D

        if SPATIAL_SCHEME == 'Finite difference':
            spatial_scheme = conservative_fd_2d
        else:
            raise NotImplementedError('Spatial discretisation not implemented.')

        if EQUATION_TYPE == 'Scalar transport':
            eqn_type = scalar_transport_2d
            eqn_flux = transport_flux
        elif EQUATION_TYPE == 'Inviscid Burgers':
            eqn_type = inviscid_burgers_2d
            eqn_flux = burgers_flux
        else:
            raise NotImplementedError('Equation type is not implemented.')

    else:
        raise NotImplementedError('Number of spatial dimensions is not supported.')


    """ Define number of grid points and CFl numbers to be tested """
    CFL_list = [0.4, 0.5, 0.8, 1.0]

    center = -0.25
    lw = 1
    ms = 3.25
    x_left = 0.
    x_right = 2. * np.pi
    y_bottom = -1.
    y_top = 1.
    fd_style = ['r.', 'gv', 'b<', 'k*']
    cab_style = ['ro', 'g^', 'b>', 'ks']
    plt.rc('text', usetex=True)

    for cfl in range(len(CFL_list)):

        """ Finite difference grid """
        var = grid(NO_OF_GRID_POINTS, DOMAIN_BOUNDARIES, BOUNDARY_CONDITIONS_TYPE)

        """ CABARET class object """
        msh = mesh([NO_OF_GRID_POINTS[0] - 1], 1, DOMAIN_BOUNDARIES, BOUNDARY_CONDITIONS_TYPE, eqn_flux)
        # Cell centre x coordinates
        x_centres = msh.x_coordinates
        # Cell faces x coordinates
        x_faces = np.linspace(DOMAIN_BOUNDARIES[0][0] - msh.cells_data[0].dx,
                              DOMAIN_BOUNDARIES[1][0] + msh.cells_data[-1].dx, NO_OF_GRID_POINTS[0] - 1 + 3)

        '''======PLEASE INPUT INITIAL CONDITIONS FOR FD======'''
        # x coordinates
        x = var.points_data.coordinates
        # Initial conditions defined with a numpy array
        ic = np.sin(x)
        # Prescribed velocity field
        v = 1. * np.ones(var.nx + 2)
        '''======================================================'''

        '''======PLEASE INPUT INITIAL CONDITIONS FOR CABARET======'''
        # Initial conditions at cell centres defined with a numpy array
        ic_centres = np.sin(x_centres)
        # Initial conditions at cell faces
        ic_faces = np.sin(x_faces)
        # Prescribed velocity field at cell centres
        v_centre = 1. * np.ones(msh.nx + 2)
        # Prescribed velocity field at cell faces
        v_faces = 1. * np.ones(msh.nx + 3)
        '''======================================================'''

        var.variables_data.values = ic
        dt = CFL_list[cfl] * var.dx / np.max(np.abs(ic))
        dt0 = CFL_list[cfl] * var.dx / np.max(np.abs(ic))

        for cell in range(msh.nx + 2):

            msh.cells_data[cell].centre[0] = ic_centres[cell]
            msh.cells_data[cell].c_centre[0] = ic_centres[cell]
            msh.cells_data[cell].dt = dt

            for j in range(msh.cells_data[cell].faces.shape[0]):
                msh.cells_data[cell].c_faces[j][0] = ic_faces[cell + j]
                msh.cells_data[cell].faces[j][0] = ic_faces[cell + j]

        msh.enforce_boundary_conditions()

        """ Plotting parameters and visualisations """
        # Initial conditions
        # plt.plot(x[1:-1], var.variables_data.values[1:-1], 'k', label='IC')
        # plt.legend(loc='best')
        # plt.ylabel(r'$\phi$')
        # # plt.axhline(H, linestyle=':', color='black')
        # plt.ylim([y_bottom, y_top])
        # plt.pause(1)

        t = 0 # Initial time
        while abs(FINAL_TIME - t) > 1.e-7:

            if t + dt >= FINAL_TIME: dt = FINAL_TIME - t

            """ Equation type """
            eqn_type(dt, v, numerical_flux, var)

            """ Spatial discretisation """
            spatial_scheme(var, 1.)

            """ Time evolution """
            temporal_scheme(var, dt)

            t += dt
            dt = CFL_list[cfl] * var.dx / np.max(np.abs(var.variables_data.values))

        dt = dt0

        t = 0  # Initial time
        while abs(FINAL_TIME - t) > 1.e-7:

            if t + dt >= FINAL_TIME: dt = FINAL_TIME - t

            """ CABARET advector """
            cabaret_1d(msh)

            t += dt
            dt = CFL_list[cfl] * msh.cells_data[0].dx / np.max(np.abs(msh.cell_centre_values[:, 0]))


        """ Numerical solution """
        plt.plot(x[1:-1], var.variables_data.values[1:-1], fd_style[cfl], label=r'FD' + str(CFL_list[cfl]), lw=lw, ms=ms, fillstyle='none')
        plt.plot(x_centres[1:-1], msh.cell_centre_values[1:-1, 0], cab_style[cfl], label=r'CAB' + str(CFL_list[cfl]), lw=lw, ms=ms, fillstyle='none')
        plt.legend(loc='lower left')


    """ Analytical solution """
    exact_solution_fd = np.sin(x - var.variables_data.values * FINAL_TIME)
    exact_solution_cab = np.sin(x_centres - msh.cell_centre_values[:, 0] * FINAL_TIME)
    plt.plot(x_centres[1:-1], exact_solution_cab[1:-1], 'k-', label='Exact', linewidth=lw)
    plt.legend(loc='best')
    plt.xlabel('$x$', usetex=True)
    plt.ylabel(r'$u$', usetex=True)
    plt.xlim([x_left, x_right])
    plt.xticks(np.linspace(0., 2. * np.pi, 5), ['0.', r'$\pi / 2$', r'$\pi$', r'$3 \pi / 2$', r'$2 \pi$'], usetex=True)
    plt.ylim([y_bottom, y_top])
    plt.title(r'Numerical vs exact solutions, ' + r'$T =$ ' + '%.2f' % FINAL_TIME + ' seconds')
    plt.grid()
    plt.tick_params(axis='both' , direction='in')
    plt.savefig(f'{path}/burgers_sine.eps', dpi=1000)
    plt.show()