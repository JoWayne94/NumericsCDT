"""
1D linear advection convergence test with smooth ICs and varying dx and CFL

Author: JWT
"""

# Import third-party libraries
import sys, os
import matplotlib.pyplot as plt
import numpy as np

# Import user settings
from conv_setup import *

# Configure system path
if (sys.platform[:3] == 'win') or (sys.platform[:3] == 'Win'):
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
else:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


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


def compute_l2_err(numerical_solution, true_solution, dx):

    return np.sqrt(np.sum(dx * np.power(numerical_solution - true_solution, 2))) / np.sqrt(np.sum(dx * np.power(true_solution, 2)))


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
    ngp_list = [101, 201, 401, 601, 801]  # testing number of cells for fair comparison with CABARET
    CFL_list = [0.4, 0.51, 0.8]

    """ Error measures """
    L2_errors_fd = np.empty((len(ngp_list), len(CFL_list)))
    L2_errors_cab = np.empty_like(L2_errors_fd)

    center = -0.25
    lw = 1.5
    ms = 5
    y_bottom = 0.
    y_top = 1.
    line_colours = ['r', 'g', 'b', 'y', 'k']

    for cfl in range(len(CFL_list)):

        for ngp in range(len(ngp_list)):

            """ Finite difference grid """
            var = grid([ngp_list[ngp]], DOMAIN_BOUNDARIES, BOUNDARY_CONDITIONS_TYPE)

            """ CABARET class object """
            msh = mesh([ngp_list[ngp] - 1], 1, DOMAIN_BOUNDARIES, BOUNDARY_CONDITIONS_TYPE, eqn_flux)
            # Cell centre x coordinates
            x_centres = msh.x_coordinates
            # Cell faces x coordinates
            x_faces = np.linspace(DOMAIN_BOUNDARIES[0][0] - msh.cells_data[0].dx,
                                  DOMAIN_BOUNDARIES[1][0] + msh.cells_data[-1].dx, ngp_list[ngp] - 1 + 3)

            '''======PLEASE INPUT INITIAL CONDITIONS FOR FD======'''
            # x coordinates
            x = var.points_data.coordinates
            # Initial conditions defined with a numpy array
            ic = np.where((x - center) % 2. < 0.5, np.power(np.sin(2. * (x - center) * np.pi), 2), 0.)
            # Prescribed velocity field
            v = 1. * np.ones(var.nx + 2)
            '''======================================================'''

            '''======PLEASE INPUT INITIAL CONDITIONS FOR CABARET======'''
            # Initial conditions at cell centres defined with a numpy array
            ic_centres = np.where((x_centres - center) % 2. < 0.5,
                                  np.power(np.sin(2. * (x_centres - center) * np.pi), 2), 0.)
            # Initial conditions at cell faces
            ic_faces = np.where((x_faces - center) % 2. < 0.5, np.power(np.sin(2. * (x_faces - center) * np.pi), 2), 0.)
            # Prescribed velocity field at cell centres
            v_centre = 1. * np.ones(msh.nx + 2)
            # Prescribed velocity field at cell faces
            v_faces = 1. * np.ones(msh.nx + 3)
            '''======================================================'''

            var.variables_data.values = ic
            dt = CFL_list[cfl] * var.dx / np.max(np.abs(v))

            for cell in range(msh.nx + 2):

                msh.cells_data[cell].centre[0] = ic_centres[cell]
                msh.cells_data[cell].c_centre[0] = v_centre[cell]
                msh.cells_data[cell].dt = dt

                for j in range(msh.cells_data[cell].faces.shape[0]):
                    msh.cells_data[cell].c_faces[j][0] = v_faces[cell + j]
                    msh.cells_data[cell].faces[j][0] = ic_faces[cell + j]  # eqn_flux(v_faces, ic_faces)[cell + j]

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

                """ CABARET advector """
                cabaret_1d(msh)

                # Replot
                # plt.cla()
                # plt.plot(x[1:-1], var.variables_data.values[1:-1], 'b', label='Time = ' + '%.2f' % (t + dt) + ' s', linewidth=lw)
                # plt.title(r'$CFL = $' + str(CFL_list[cfl]))
                # plt.legend(loc='lower left')
                # plt.xlabel('x')
                # plt.ylabel(r'$\phi$')
                # plt.ylim([y_bottom, y_top])
                # plt.pause(0.01)

                t += dt

            """ Analytical solution """
            exact_solution_fd = np.where((x[1:-1] - center) % 2. < 0.5, np.power(np.sin(2. * (x[1:-1] - center) * np.pi), 2), 0.)
            exact_solution_cab = np.where((x_centres[1:-1] - center) % 2. < 0.5, np.power(np.sin(2. * (x_centres[1:-1] - center) * np.pi), 2), 0.)
            # plt.plot(x[1:-1], exact_solution, 'k--', label='Exact soln', linewidth=lw)
            # plt.legend(loc='best')
            # plt.ylabel(r'$\phi$')
            # plt.title(r'$CFL = $' + str(CFL) + ', ' + r'$T = L^{-1} u$')
            # plt.axhline(0, linestyle=':', color='black')
            # plt.grid()
            # plt.show()

            """ Append L2 error to list """
            L2_errors_fd[ngp][cfl] = compute_l2_err(var.variables_data.values[1:-1], exact_solution_fd, var.dx)
            L2_errors_cab[ngp][cfl] = compute_l2_err(msh.cell_centre_values[1:-1, 0], exact_solution_cab, var.dx)

    # print(L2_errors)

    """ Plot convergence """
    for c in range(len(CFL_list)):

        plt.plot(np.log10(np.array(ngp_list) - 1), np.log10(L2_errors_fd[:, c]), 'x--',
                 label='1stOrder, CFL = ' + str(CFL_list[c]), linewidth=lw, color=line_colours[c])

        plt.plot(np.log10(np.array(ngp_list) - 1), np.log10(L2_errors_cab[:, c]), 'o--',
                 label='2ndOrder, CFL = ' + str(CFL_list[c]), linewidth=lw, markersize=ms, fillstyle='none', color=line_colours[c])

    # plt.rc('text', usetex=True)
    # plt.rc('font', family='serif')

    # Define the base of the triangle on the log-log scale
    base_x = 2.35  # Starting x-coordinate
    base_y = -0.8  # Starting y-coordinate
    length = 0.2  # Length of the triangle sides on log scale
    gradient = -1  # Desired slope of the triangle on log-log scale

    # Calculate the triangle vertices for the general gradient
    # Top vertex: move `length` along x and y axes to maintain slope of 1
    x_vertices = [base_x, base_x + length, base_x + length]
    y_vertices = [base_y, base_y, base_y + gradient * length]

    # Plot the triangle
    plt.plot(x_vertices + [x_vertices[0]], y_vertices + [y_vertices[0]], 'k', linewidth=lw)
    plt.text(2.475, -0.9, r'\textbf{m = 1}', usetex=True)

    # Define the base of the triangle on the log-log scale
    base_x = 2.3  # Starting x-coordinate
    base_y = -2.3  # Starting y-coordinate
    length = 0.2  # Length of the triangle sides on log scale
    gradient = -2  # Desired slope of the triangle on log-log scale

    # Calculate the triangle vertices for the general gradient
    # Top vertex: move `length` along x and y axes to maintain slope of 1
    x_vertices = [base_x, base_x + length, base_x + length]
    y_vertices = [base_y, base_y, base_y + gradient * length]

    # Plot the triangle
    plt.plot(x_vertices + [x_vertices[0]], y_vertices + [y_vertices[0]], 'k', linewidth=lw)
    plt.text(2.425, -2.5, r'\textbf{m = 2}', usetex=True)

    plt.legend(loc='best')
    plt.xlabel(r'$n_x$')
    plt.xticks(np.log10(np.array(ngp_list) - 1), [str(n - 1) for n in ngp_list])
    plt.ylabel(r'$L_2$ error')
    plt.title(r'Convergence on a $log_{10}-log_{10}$ scale')  # , $T = $' + str(FINAL_TIME) + ' seconds'
    # plt.axhline(0, linestyle=':', color='black')
    plt.grid()
    plt.tick_params(axis='both' , direction='in')
    plt.show()
