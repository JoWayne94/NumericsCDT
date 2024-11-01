"""
A top-hat translating to the left or right

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

# Import CABARET methods
from src.CABARET.CABARET import *
from src.CABARET.mesh.mesh1D import *

# Import equation types
from src.solvers.scalartransport import *
from src.solvers.inviscidburgers import *


if __name__ == '__main__':
    """
    main()
    """

    if EQUATION_TYPE == 'Scalar transport':
        eqn_flux = transport_flux
    elif EQUATION_TYPE == 'Inviscid Burgers':
        eqn_flux = burgers_flux
    else:
        raise NotImplementedError('Equation type is not implemented.')

    if NO_OF_DIMENSIONS == 1:
        mesh = Mesh1D
    elif NO_OF_DIMENSIONS == 2:
        mesh = Mesh1D
    else:
        raise NotImplementedError('Number of spatial dimensions is not supported.')

    """ Mesh class object """
    msh = mesh(NO_OF_CELLS, NO_OF_VARIABLES, DOMAIN_BOUNDARIES, BOUNDARY_CONDITIONS_TYPE, eqn_flux)
    # Cell centre x coordinates
    x_centres = msh.x_coordinates
    # Cell faces x coordinates
    x_faces = np.linspace(DOMAIN_BOUNDARIES[0][0] - msh.cells_data[0].dx,
                          DOMAIN_BOUNDARIES[1][0] + msh.cells_data[-1].dx, NO_OF_CELLS[0] + 3)

    '''======PLEASE INPUT YOUR INITIAL CONDITIONS BELOW======'''
    # Initial conditions at cell centres defined with a numpy array
    ic_centres = np.where((0.5 <= x_centres) & (x_centres <= 1.5), 1., 0.)
    # Initial conditions at cell faces
    ic_faces = np.where((0.5 <= x_faces) & (x_faces <= 1.5), 1., 0.)
    # Prescribed velocity field at cell centres
    v_centre = 1. * np.ones(msh.nx + 2)
    # Prescribed velocity field at cell faces
    v_faces = 1. * np.ones(msh.nx + 3)
    '''======================================================'''

    # For asynchronous time-stepping in the future, dt will be local
    dt = CFL * msh.cells_data[0].dx / np.max(np.abs(v_centre))

    for cell in range(msh.nx + 2):

        msh.cells_data[cell].centre[0] = ic_centres[cell]
        msh.cells_data[cell].c_centre[0] = v_centre[cell]
        msh.cells_data[cell].dt = dt

        for j in range(msh.cells_data[cell].faces.shape[0]):
            msh.cells_data[cell].c_faces[j][0] = v_faces[cell + j]
            msh.cells_data[cell].faces[j][0] = ic_faces[cell + j]  # eqn_flux(v_faces, ic_faces)[cell + j]

    # Python objects are equated with pointers, so BCs only need to be initialised once at the beginning
    msh.enforce_boundary_conditions()

    """ Plotting parameters and visualisations """
    lw = 1.5
    y_bottom = -0.5
    y_top = 1.5

    # Initial conditions
    plt.plot(x_centres[1:-1], msh.cell_centre_values[1:-1, 0], 'k', label='IC')
    plt.legend(loc='best')
    plt.ylabel(r'$\phi$')
    # plt.axhline(H, linestyle=':', color='black')
    plt.ylim([y_bottom, y_top])
    plt.pause(1)

    t = 0  # Initial time
    while abs(FINAL_TIME - t) > 1.e-7:

        if t + dt >= FINAL_TIME: dt = FINAL_TIME - t

        """ CABARET advector """
        cabaret_1d(msh)

        # Replot
        plt.cla()
        plt.plot(x_centres[1:-1], msh.cell_centre_values[1:-1, 0], 'b', label='Time = ' + '%.2f' % (t + dt) + ' s',
                 linewidth=lw)
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
    exact_solution = np.where((0.5 <= x_centres[1:-1]) & (x_centres[1:-1] <= 1.5), 1., 0.)
    plt.plot(x_centres[1:-1], exact_solution, 'k--', label='Exact soln', linewidth=lw)
    plt.legend(loc='best')
    plt.ylabel(r'$\phi$')
    plt.title(r'$CFL = $' + str(CFL) + ', ' + r'$T = L^{-1} u$')
    plt.axhline(0, linestyle=':', color='black')
    plt.grid()
    plt.show()
