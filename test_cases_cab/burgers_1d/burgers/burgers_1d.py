"""
Inviscid Burgers' equation test case

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
    ic_centres = np.empty_like(x_centres)
    # Initial conditions at cell faces
    ic_faces = np.empty_like(x_faces)

    for i in range(x_centres.shape[0]):
        if x_centres[i] < 0.: ic_centres[i] = 0.
        if 0. <= x_centres[i] < 1.: ic_centres[i] = 1.
        else: ic_centres[i] = 0.

    for i in range(x_faces.shape[0]):
        if x_faces[i] < 0.: ic_faces[i] = 0.
        if 0. <= x_faces[i] < 1.: ic_faces[i] = 1.
        else: ic_faces[i] = 0.
    '''======================================================'''

    dt = CFL * msh.cells_data[0].dx / np.max(np.abs(ic_centres))

    for cell in range(msh.nx + 2):

        msh.cells_data[cell].centre[0] = ic_centres[cell]
        msh.cells_data[cell].c_centre[0] = ic_centres[cell]
        msh.cells_data[cell].dt = dt

        for j in range(msh.cells_data[cell].faces.shape[0]):
            msh.cells_data[cell].c_faces[j][0] = ic_faces[cell + j]
            msh.cells_data[cell].faces[j][0] = ic_faces[cell + j]

    msh.enforce_boundary_conditions()

    """ Plotting parameters and visualisations """
    lw = 1.5
    y_bottom = -0.2
    y_top = 1.5

    # Initial conditions
    plt.plot(x_centres[1:-1], msh.cell_centre_values[1:-1, 0], 'k', label='IC')
    plt.legend(loc='best')
    plt.ylabel(r'$u$')
    # plt.axhline(H, linestyle=':', color='black')
    plt.ylim([y_bottom, y_top])
    plt.pause(1)

    t = 0 # Initial time
    while abs(FINAL_TIME - t) > 1.e-7:

        if t + dt >= FINAL_TIME: dt = FINAL_TIME - t

        """ CABARET advector """
        cabaret_1d(msh)

        # Replot
        plt.cla()
        plt.plot(x_centres[1:-1], msh.cell_centre_values[1:-1, 0], 'bo', label='Time = ' + '%.2f' % (t + dt) + ' s', linewidth=lw)
        plt.title(r'$CFL = $' + str(CFL))
        plt.legend(loc='lower left')
        plt.xlabel('x')
        plt.ylabel(r'$u$')
        plt.ylim([y_bottom, y_top])
        plt.pause(0.01)

        t += dt
        dt = CFL * msh.cells_data[0].dx / np.max(np.abs(msh.cell_centre_values[:, 0]))


    """ Analytical solution """
    ac = np.zeros_like(x_centres)
    for i in range(x_centres.shape[0]):
        if 0. <= x_centres[i] < FINAL_TIME:
            ac[i] = float(x_centres[i]) / FINAL_TIME
        if FINAL_TIME <= x_centres[i] < 0.5 * FINAL_TIME + 1.:
            ac[i] = 1.

    plt.plot(x_centres[1:-1], ac[1:-1], 'k--', label='Exact soln', linewidth=lw)
    plt.legend(loc='best')
    plt.ylabel(r'$u$')
    plt.title(r'$CFL = $' + str(CFL) + ', ' + r'$T = $' + str(FINAL_TIME) + r'$ s$')
    plt.axhline(0, linestyle=':', color='black')
    plt.grid()
    plt.show()
