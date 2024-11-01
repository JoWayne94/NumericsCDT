# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
"""
Main module

Note:

    1. Serves as a template to create other solvers
    2. A second-order CABARET is implemented
    3. Computational domain only limited to quadrilaterals
    4. Grid spacings are constant for each spatial dimension (\delta x is not a function of x)

Future works:

1.

Author: JWT
"""

# Import third-party libraries
import sys, os
import matplotlib.pyplot as plt

# Import user settings
from setup_cab import *

# Import CABARET methods
from src.CABARET.CABARET import *

# Import grid classes
from src.CABARET.mesh.mesh1D import *

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
    x_faces = np.linspace(DOMAIN_BOUNDARIES[0][0] - msh.cells_data[0].dx, DOMAIN_BOUNDARIES[1][0] + msh.cells_data[-1].dx, NO_OF_CELLS[0] + 3)

    '''======PLEASE INPUT YOUR INITIAL CONDITIONS BELOW======'''
    center = -0.25
    # Initial conditions at cell centres defined with a numpy array
    ic_centres = np.where((x_centres - center) % 2. < 0.5, np.power(np.sin(2. * (x_centres - center) * np.pi), 2), 0.)
    # Initial conditions at cell faces
    ic_faces = np.where((x_faces - center) % 2. < 0.5, np.power(np.sin(2. * (x_faces - center) * np.pi), 2), 0.)
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
    y_bottom = 0.
    y_top = 1.

    # Initial conditions
    plt.plot(x_centres[1:-1], msh.cell_centre_values[1:-1, 0], 'k', label='IC')
    plt.legend(loc='best')
    plt.ylabel(r'$\phi$')
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
        plt.plot(x_centres[1:-1], msh.cell_centre_values[1:-1, 0], 'b', label='Time = ' + '%.2f' % (t + dt) + ' s', linewidth=lw)
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
    exact_solution = np.where((x_centres[1:-1] - center) % 2. < 0.5, np.power(np.sin(2. * (x_centres[1:-1] - center) * np.pi), 2), 0.)
    plt.plot(x_centres[1:-1], exact_solution, 'k--', label='Exact soln', linewidth=lw)
    plt.legend(loc='best')
    plt.ylabel(r'$\phi$')
    plt.title(r'$CFL = $' + str(CFL) + ', ' + r'$T = L u^{-1}$')
    plt.axhline(0, linestyle=':', color='black')
    plt.grid()
    plt.show()
