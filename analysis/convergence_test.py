"""
1D linear advection convergence test with smooth ICs (sine wave) and varying dx and CFL

Note: Takes about 3 minutes to run

Author: JWT
"""

# Import third-party libraries
import sys, os
import matplotlib.pyplot as plt
import numpy as np
from itertools import combinations

# Import user settings
from conv_setup import *

# Configure system path
path = os.path.dirname(__file__)

if (sys.platform[:3] == 'win') or (sys.platform[:3] == 'Win'):
    sys.path.append(os.path.abspath(os.path.join(path, '..')))
else:
    sys.path.append(os.path.abspath(os.path.join(path, '..')))


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
    # ngp_list = [65, 129, 257, 513, 1025]  # testing number of cells for fair comparison with CABARET
    # CFL_list = [0.4, 0.5, 0.8]
    ngp_list = []
    CFL_list = []

    """ Error measures """
    L2_errors_fd = np.empty((len(ngp_list), len(CFL_list)))
    L2_errors_cab = np.empty_like(L2_errors_fd)

    center = -0.25
    lw = 1.5
    ms = 5
    y_bottom = 0.
    y_top = 1.
    line_colours = ['rx--', 'go--', 'bs--', 'y', 'k']

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
            dt = CFL_list[cfl] * var.dx / np.max(np.abs(v))

            for cell in range(msh.nx + 2):

                msh.cells_data[cell].centre[0] = ic_centres[cell]
                msh.cells_data[cell].c_centre[0] = v_centre[cell]
                msh.cells_data[cell].dt = dt

                for j in range(msh.cells_data[cell].faces.shape[0]):
                    msh.cells_data[cell].c_faces[j][0] = v_faces[cell + j]
                    msh.cells_data[cell].faces[j][0] = ic_faces[cell + j]  # eqn_flux(v_faces, ic_faces)[cell + j]

            msh.enforce_boundary_conditions()

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

                t += dt

            """ Analytical solution """
            exact_solution_fd = np.sin(x[1:-1] - v[1:-1] * FINAL_TIME)
            exact_solution_cab = np.sin(x_centres[1:-1] - v_centre[1:-1] * FINAL_TIME)

            """ Append L2 error to list """
            L2_errors_fd[ngp][cfl] = compute_l2_err(var.variables_data.values[1:-1], exact_solution_fd, var.dx)
            L2_errors_cab[ngp][cfl] = compute_l2_err(msh.cell_centre_values[1:-1, 0], exact_solution_cab, var.dx)


    """ Manual input of errors for faster plotting """
    ngp_list = [65., 129., 257., 513., 1025.]
    CFL_list = [0.4, 0.5, 0.8]
    dx_list = [(DOMAIN_BOUNDARIES[1][0] - DOMAIN_BOUNDARIES[0][0]) / d for d in np.array(ngp_list) - 1]

    L2_errors_fd = np.array([[0.16899958, 0.1429633, 0.05983103],
                            [0.08838383, 0.07421572, 0.03037298],
                            [0.04521095, 0.03782036, 0.01530311],
                            [0.02286653, 0.01909208, 0.007681],
                            [0.01149933, 0.009592, 0.00384789]])

    L2_errors_cab = np.array([[1.85288587e-03, 6.41224650e-08, 1.26143326e-03],
                              [5.48948872e-04, 2.83469309e-09, 3.41181834e-04],
                              [1.69172889e-04, 1.25286357e-10, 9.08867718e-05],
                              [5.17634909e-05, 5.53722377e-12, 2.39328136e-05],
                              [1.60151578e-05, 2.44871448e-13, 6.18840552e-06]])

    """ Compute numerical order of convergence """
    # for cfl in range(len(CFL_list)):
    #
    #     # Calculate n for each unique pair of (epsilon, delta x) values
    #     n_values = []
    #     for (i, j) in combinations(range(len(L2_errors_fd[:, cfl])), 2):
    #         epsilon_1, epsilon_2 = L2_errors_fd[:, cfl][i], L2_errors_fd[:, cfl][j]
    #         delta_x_1, delta_x_2 = dx_list[i], dx_list[j]
    #
    #         # Using the formula provided to calculate n for each pair
    #         n = (np.log(epsilon_1) - np.log(epsilon_2)) / (np.log(delta_x_1) - np.log(delta_x_2))
    #         n_values.append(n)
    #
    #     # Calculate the average n
    #     average_n = np.mean(n_values)
    #
    #     # Output the result
    #     print(f"Calculated n values for each pair for CFL = {CFL_list[cfl]}:", n_values)
    #     print(f"Average value of n for CFL = {CFL_list[cfl]}:", average_n)

    plt.rc('text', usetex=True)
    """ Plot convergence for FD """
    plt.xscale('log')
    plt.yscale('log')

    for c in range(len(CFL_list)):

        plt.plot(np.array(ngp_list) - 1, L2_errors_fd[:, c], line_colours[c],
                 label='FD' + str(CFL_list[c]), lw=lw, ms=ms, fillstyle='none')

    # Define the desired gradient and constant
    m1 = -1  # Desired gradient (slope) on the log-log scale
    m2 = -2
    C1 = 7  # Constant that sets the vertical position of the line
    C2 = 3

    # Calculate y based on the power law y = C * x^m
    y1 = C1 * (np.array(ngp_list) - 1) ** m1
    y2 = C2 * (np.array(ngp_list) - 1) ** m2

    plt.plot(np.array(ngp_list) - 1, y1, 'k:', label=r'$A \Delta x$', lw=lw, ms=ms)

    plt.legend(loc='best')
    plt.xlabel(r'$\Delta x$', usetex=True)
    plt.xticks(np.array(ngp_list) - 1, [f"$2 \pi / {int(nx)}$" for nx in np.array(ngp_list) - 1])
    plt.ylabel(r'$L_2$ error', usetex=True)
    plt.title(r'Convergence on a log$_{10}$-log$_{10}$ scale', usetex=True)  # , $T = $' + str(FINAL_TIME) + ' seconds'
    plt.grid()
    # plt.tick_params(axis='both', direction='in')
    plt.savefig(f'{path}/convergence_fd.eps', dpi=1000)
    plt.show()

    """ Plot convergence for CABARET """
    # for cfl in range(len(CFL_list)):
    #
    #     # Calculate n for each unique pair of (epsilon, delta x) values
    #     n_values = []
    #     for (i, j) in combinations(range(len(L2_errors_cab[:, cfl])), 2):
    #         epsilon_1, epsilon_2 = L2_errors_cab[:, cfl][i], L2_errors_cab[:, cfl][j]
    #         delta_x_1, delta_x_2 = dx_list[i], dx_list[j]
    #
    #         # Using the formula provided to calculate n for each pair
    #         n = (np.log(epsilon_1) - np.log(epsilon_2)) / (np.log(delta_x_1) - np.log(delta_x_2))
    #         n_values.append(n)
    #
    #     # Calculate the average n
    #     average_n = np.mean(n_values)
    #
    #     # Output the result
    #     print(f"Calculated n values for each pair for CFL = {CFL_list[cfl]}:", n_values)
    #     print(f"Average value of n for CFL = {CFL_list[cfl]}:", average_n)

    plt.clf()
    plt.xscale('log')
    plt.yscale('log')

    for c in range(len(CFL_list)):

        plt.plot(np.array(ngp_list) - 1, L2_errors_cab[:, c], line_colours[c],
                 label='CAB' + str(CFL_list[c]), linewidth=lw, markersize=ms, fillstyle='none')

    plt.plot(np.array(ngp_list) - 1, y2, 'k:', lw=lw, ms=ms, label=r'$A \Delta x^{2}$')

    plt.legend(loc='best')
    plt.xlabel(r'$\Delta x$', usetex=True)
    plt.xticks(np.array(ngp_list) - 1, [f"$2 \pi / {int(nx)}$" for nx in np.array(ngp_list) - 1])
    plt.ylabel(r'$L_2$ error', usetex=True)
    plt.title(r'Convergence on a log$_{10}$-log$_{10}$ scale', usetex=True)
    plt.grid()
    # plt.tick_params(axis='both', direction='in')
    plt.savefig(f'{path}/convergence_cab.eps', dpi=1000)
    plt.show()
